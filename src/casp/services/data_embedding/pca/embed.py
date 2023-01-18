import argparse
import pickle
import time

import pandas as pd
from google.cloud import storage
from torch.utils.data import DataLoader

from casp.ml.data.dataset import CASDataset
from casp.ml.dump_manager import DumpManager
from casp.utils import get_google_service_credentials


def get_dump_manager(bucket, filepath) -> "DumpManager":
    blob = bucket.blob(filepath)
    pickle_in = blob.download_as_string()
    return pickle.loads(pickle_in)


def write_embeddings(bucket, embeddings, ids, output_path, postfix):
    df = pd.DataFrame(embeddings.cpu().numpy())
    cas_id_values = ids.cpu().numpy().astype(int)
    df.insert(0, "db_ids", cas_id_values)
    bucket.blob(f"{output_path}/batch_{postfix}.csv").upload_from_string(
        df.to_csv(header=False, index=False), "text/csv"
    )


def main(bucket_name, data_storage_path, dm_storage_path, output_storage_path):
    credentials, project_id = get_google_service_credentials()
    storage_client = storage.Client(project=project_id, credentials=credentials)
    bucket = storage_client.bucket(bucket_name=bucket_name)
    dm = get_dump_manager(bucket, dm_storage_path)
    dm.to_cuda()

    model = dm.model
    one_pass = dm.running_stat
    t = dm.transform

    dataset = CASDataset(
        storage_client=storage_client,
        storage_path=data_storage_path,
        bucket_name=bucket_name,
        use_gpu=True,
        transform=t,
        running_stat=one_pass,
    )
    dataloader = DataLoader(dataset=dataset, batch_size=10000)

    start = time.time()
    counter = 1
    print("===== START EMBEDDING =====")
    for batch, ids in dataloader:
        embeddings = model.transform(batch)
        write_embeddings(
            bucket=bucket, embeddings=embeddings, ids=ids, output_path=output_storage_path, postfix=counter
        )
        if counter % 25 == 0:
            print(f"----- Processed {counter} chunks {time.time() - start}")
            start = time.time()

        counter += 1

    print("===== FINISH EMBEDDING =====")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--bucket_name",
        help="Bucket name to work with. This bucket will be use to extract data and for saving a model",
    )
    parser.add_argument("--data_storage_path", help="Location in the bucket where to extract the data from")
    parser.add_argument("--dm_storage_path", help="Location in the bucket where to get a DumpManager pickle")
    parser.add_argument("--output_storage_path", help="Location in the bucket where to write output embeddings")
    args = parser.parse_args()

    main(
        bucket_name=args.bucket_name,
        data_storage_path=args.data_storage_path,
        dm_storage_path=args.dm_storage_path,
        output_storage_path=args.output_storage_path,
    )

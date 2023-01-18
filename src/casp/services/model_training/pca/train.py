import argparse
import time

from google.cloud import storage
from torch.utils.data import DataLoader

from casp.ml.data import transforms
from casp.ml.data.dataset import CASDataset
from casp.ml.dump_manager import DumpManager
from casp.ml.models.pca import IncrementalPCA
from casp.ml.running_stats import OnePassMeanVarStd
from casp.services.utils import get_google_service_credentials


def save_check_point(model, running_stat, t, bucket, postfix):
    dm = DumpManager(model=model, running_stat=running_stat, transform=t)
    dm.to_cpu()
    f_name = f"./dump_manager_{postfix}.pickle"
    dm.save(f_name)
    dm_blob = bucket.blob(f"models/dump_manager_{postfix}.pickle")
    dm_blob.upload_from_filename(f_name)
    dm.to_cuda()


def main(bucket_name, storage_path):
    credentials, project_id = get_google_service_credentials()
    storage_client = storage.Client(project=project_id, credentials=credentials)
    bucket = storage_client.bucket(bucket_name=bucket_name)

    model = IncrementalPCA(n_components=512)
    one_pass = OnePassMeanVarStd()
    t = transforms.Compose(
        [
            transforms.RowWiseNormalization(sc_rna_normalization_pseudocount=10000),
            transforms.ColumnWiseNormalization(one_pass_mean_var_std=one_pass),
        ]
    )
    dataset = CASDataset(
        storage_client=storage_client,
        storage_path=storage_path,
        bucket_name=bucket_name,
        use_gpu=True,
        transform=t,
        running_stat=one_pass,
    )
    dataset.calculate_running_stats()
    dataloader = DataLoader(dataset=dataset, batch_size=5000)

    start = time.time()
    counter = 1
    print("===== START TRAINING =====")
    for batch, ids in dataloader:
        model(batch)
        if counter % 25 == 0:
            save_check_point(model=model, running_stat=one_pass, t=t, bucket=bucket, postfix=counter)
            print(f"----- Processed {counter} chunks {time.time() - start}")
            start = time.time()

        counter += 1

    save_check_point(model=model, running_stat=one_pass, t=t, bucket=bucket, postfix="final")
    print("===== FINISH TRAINING =====")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--bucket_name",
        help="Bucket name to work with. This bucket will be use to extract data and for saving a model",
    )
    parser.add_argument("--data_storage_path", help="Location in the bucket where to extract the data from")
    args = parser.parse_args()
    main(bucket_name=args.bucket_name, storage_path=args.data_storage_path)

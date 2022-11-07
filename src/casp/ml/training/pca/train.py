from pathlib import Path
import time
import torch
from torch.utils.data import DataLoader
from google.cloud import storage
from casp.ml.training.pca.utils import get_google_service_credentials
from casp.ml.models.pca import IncrementalPCA
from casp.ml.running_stats import OnePassMeanVarStd
from casp.ml.data import transforms
from casp.ml.data.dataset import CASDataset
from casp.ml.dump_manager import DumpManager


def save_check_point(model, running_stat, t, bucket, counter):
    dm = DumpManager(model=model, running_stat=running_stat, transform=t)
    dm.to_cpu()
    f_name = f"./dump_manager_{counter}.pickle"
    dm.save(f_name)
    dm_blob = bucket.blob(f"models/dump_manager_{counter}.pickle")
    dm_blob.upload_from_filename(f_name)
    dm.to_cuda()


def main():
    credentials, project_id = get_google_service_credentials()
    storage_client = storage.Client(project=project_id, credentials=credentials)
    bucket_name = "fedor-test-bucket"
    storage_path = "test_data"
    bucket = storage_client.bucket(bucket_name=bucket_name)

    model = IncrementalPCA(n_components=512)
    one_pass = OnePassMeanVarStd()
    t = transforms.Compose(
        [
            transforms.RowWiseNormalization(sc_rna_normalization_pseudocount=10000),
            transforms.ColumnWiseNormalization(one_pass_mean_var_std=one_pass)
        ]
    )
    dataset = CASDataset(
        storage_client=storage_client,
        storage_path=storage_path,
        bucket_name=bucket_name,
        use_gpu=True,
        transform=t,
        running_stat=one_pass
    )
    dataloader = DataLoader(dataset=dataset, batch_size=10000)

    start = time.time()
    counter = 1
    print("===== START TRAINING =====")
    for batch, ids in dataloader:
        model(batch)
        start = time.time()

        if counter % 25 == 0:
            save_check_point(model=model, running_stat=one_pass, t=t, bucket=bucket, counter=counter)
            print(f"Processed {counter} chunks", time.time() - start, sep="\n")
            start = time.time()

        counter += 1


if __name__ == "__main__":
    print("I AM STARTING!")
    Path("./data").mkdir(parents=True, exist_ok=True)
    main()

import argparse
import time

import neptune
from google.cloud import storage
from torch.utils.data import DataLoader

from casp.services import settings
from casp.ml.data import transforms
from casp.ml.data.dataset import CASDataset
from casp.ml.dump_manager import DumpManager
from casp.ml.models.pca import IncrementalPCA
from casp.ml.running_stats import OnePassMeanVarStd
from casp.services.utils import get_google_service_credentials


def save_check_point(model, running_stat, t, bucket, postfix, save_path, use_gpu):
    dm = DumpManager(model=model, running_stat=running_stat, transform=t)
    dm.to_cpu()
    f_name = f"./dump_manager_{postfix}.pickle"
    dm.save(f_name)
    dm_blob = bucket.blob(f"{save_path}/dump_manager_{postfix}.pickle")
    dm_blob.upload_from_filename(f_name)
    if use_gpu:
        dm.to_cuda()


def main(bucket_name, storage_path, checkpoint_save_path: str, use_gpu: bool, batch_size: int, n_components: int):
    credentials, project_id = get_google_service_credentials()
    experiment = neptune.init_run(
        project="fedorgrab/cas",
        name="Incremental PCA 4m cells",
        api_token=settings.NEPTUNE_API_KEY,
    )
    storage_client = storage.Client(project=project_id, credentials=credentials)
    bucket = storage_client.bucket(bucket_name=bucket_name)

    model = IncrementalPCA(n_components=n_components)
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
        use_gpu=use_gpu,
        transform=t,
        running_stat=one_pass,
    )
    dataset.calculate_running_stats()
    dataloader = DataLoader(dataset=dataset, batch_size=batch_size)

    start = time.time()
    counter = 1
    print("===== START TRAINING =====")
    for batch, ids in dataloader:
        model(batch)
        if counter % 100 == 0:
            save_check_point(
                model=model, 
                running_stat=one_pass,
                t=t,
                bucket=bucket,
                postfix=counter,
                save_path=checkpoint_save_path, 
                use_gpu=use_gpu
            )
            print(f"----- Processed {counter} chunks {time.time() - start}")
            start = time.time()

        counter += 1

    save_check_point(
        model=model,
        running_stat=one_pass,
        t=t, bucket=bucket,
        postfix="final",
        save_path=checkpoint_save_path,
        use_gpu=use_gpu
    )
    print("===== FINISH TRAINING =====")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--bucket_name",
        help="Bucket name to work with. This bucket will be use to extract data and for saving a model",
    )
    parser.add_argument("--data_storage_path", help="Location in the bucket where to extract the data from")
    parser.add_argument("--checkpoint_save_path", help="Directory where to save a model")
    parser.add_argument("--use_gpu", type=bool, default=False, required=False)
    parser.add_argument("--batch_size", type=int)
    parser.add_argument("--n_components", type=int)
    args = parser.parse_args()

    main(
        bucket_name=args.bucket_name,
        storage_path=args.data_storage_path,
        checkpoint_save_path=args.checkpoint_save_path,
        use_gpu=args.use_gpu,
        batch_size=args.batch_size,
        n_components=args.n_components,
    )

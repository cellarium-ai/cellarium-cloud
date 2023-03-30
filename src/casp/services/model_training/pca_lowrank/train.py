import argparse
import time

import matplotlib.pyplot as plt
import neptune
import numpy as np
import pandas as pd
import torch
# import umap
from google.cloud import bigquery, storage
from sklearn.preprocessing import StandardScaler
from sklearn.manifold import TSNE
from torch.utils.data import DataLoader

from casp.ml.data import transforms
from casp.ml.data.dataset import CASDataset
from casp.ml.dump_manager import DumpManager
from casp.ml.models.pca_lowrank import LowRankIncrementalPCA
from casp.ml.running_stats import OnePassMeanVarStd
from casp.services import settings
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


def get_celltypes(cas_ids):
    credentials, project_id = get_google_service_credentials()
    bq_client = bigquery.Client(project=project_id, credentials=credentials)
    cas_ids = [str(x) for x in cas_ids.int().tolist()]
    query = f"""
    select cell_type from cas_benchmark_v1.cas_cell_info
                 where cas_cell_index in ({", ".join(cas_ids)});
    """
    results = bq_client.query(query).result()
    cell_types = [x.cell_type for x in results]
    cell_types_simplified = []
    for ct in cell_types:
        if "T cell" in ct:
            cell_types_simplified.append("T cell")
        elif "thymocyte" in ct:
            cell_types_simplified.append("thymocyte")
        elif "B cell" in ct:
            cell_types_simplified.append("B cell")
        elif "Kupffer cell" in ct:
            cell_types_simplified.append("Kuppfer cell")
        elif "Mueler cell" in ct:
            cell_types_simplified.append("Mueler cell")
        elif "acinar cell" in ct:
            cell_types_simplified.append("acinar cell")
        elif "basal cell" in ct:
            cell_types_simplified.append("basal cell")
        elif "muscle cell" in ct:
            cell_types_simplified.append("muscle cell")
        elif "erythrocyte" in ct:
            cell_types_simplified.append("erythrocyte")
        elif "fibroblast" in ct:
            cell_types_simplified.append("fibroblast")
        elif "natural killer" in ct or "NK" in ct:
            cell_types_simplified.append("natural killer")
        elif "endothelial" in ct:
            cell_types_simplified.append("endothelial cell")
        elif "epithelial" in ct:
            cell_types_simplified.append("epithelial cell")
        elif "macrophage" in ct:
            cell_types_simplified.append("macrophage")
        else:
            cell_types_simplified.append("other")

    return cell_types_simplified


def train(dataloader, model, one_pass, t, bucket, checkpoint_save_path, use_gpu):
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
                use_gpu=use_gpu,
            )
            print(f"----- Processed {counter} chunks {time.time() - start}")
            start = time.time()

        counter += 1


def evaluation(dataloader, model, experiment):
    embeddings = []
    cell_types = []

    c = 1
    for batch, cas_ids in dataloader:
        embeddings.append(model.transform(batch))
        cell_types.extend(get_celltypes(cas_ids))

        if c == 2:
            break
        c += 1

    embeddings = torch.vstack(embeddings).cpu()
    cell_types = np.array(cell_types)
    ct_unique = np.unique(cell_types)
    # reducer = umap.UMAP()
    reducer = TSNE()
    # scaled_embeddings = StandardScaler().fit_transform(embeddings)
    scaled_embeddings = embeddings
    embedding = reducer.fit_transform(scaled_embeddings)
    figure, ax = plt.subplots(1, 1, figsize=(5, 5))

    for i, ct in enumerate(ct_unique):
        ax.scatter(
            x=embedding[cell_types == ct, 0],
            y=embedding[cell_types == ct, 1],
            label=ct_unique[i],
            s=1,
        )
    ax.legend()

    t_np = embeddings.numpy()
    df = pd.DataFrame(t_np)
    df.to_csv("embeddings_eval.tsv", index=False, sep="\t")

    df = pd.DataFrame(cell_types)
    df.to_csv("embeddings_eval_meta.tsv", index=False, sep="\t")

    experiment["embeddings_eval"].upload("embeddings_eval.tsv")
    experiment["embeddings_eval_meta"].upload("embeddings_eval_meta.tsv")
    # experiment["umap-plot"].upload(figure)
    experiment["tsne-plot"].upload(figure)


def main(bucket_name, storage_path, checkpoint_save_path, n_components, q, batch_size, use_gpu=False):
    credentials, project_id = get_google_service_credentials()
    storage_client = storage.Client(project=project_id, credentials=credentials)
    bucket = storage_client.bucket(bucket_name=bucket_name)

    experiment = neptune.init_run(
        project="fedorgrab/cas",
        api_token=settings.NEPTUNE_API_KEY,
    )
    model = LowRankIncrementalPCA(n_components=n_components, q=q)
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
    dataloader = DataLoader(dataset=dataset, batch_size=batch_size)
    train(dataloader, model, one_pass, t, bucket, checkpoint_save_path, use_gpu)
    save_check_point(
        model=model,
        running_stat=one_pass,
        t=t,
        bucket=bucket,
        postfix="final",
        save_path=checkpoint_save_path,
        use_gpu=use_gpu,
    )
    print("===== FINISH TRAINING =====")
    evaluation(dataloader, model, experiment)
    print("===== FINISH EVALUATION =====")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--bucket_name",
        help="Bucket name to work with. This bucket will be use to extract data and for saving a model",
    )
    parser.add_argument("--data_storage_path", help="Location in the bucket where to extract the data from")
    parser.add_argument("--checkpoint_save_path", help="Location in the bucket where to save the trained model")
    parser.add_argument("--n_components", help="Number of PC components", type=int)
    parser.add_argument("--q", help="Number of features for approximation matrix", type=int)
    parser.add_argument("--batch_size", help="Batch size of each iteration", type=int)
    parser.add_argument("--use_gpu", default=False, help="A boolean flag pointing whether or not ot use GPU")
    args = parser.parse_args()
    main(
        bucket_name=args.bucket_name,
        storage_path=args.data_storage_path,
        checkpoint_save_path=args.checkpoint_save_path,
        batch_size=args.batch_size,
        n_components=args.n_components,
        q=args.q,
        use_gpu=args.use_gpu,
    )

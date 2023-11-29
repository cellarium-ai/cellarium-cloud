import io
import multiprocessing
import typing as t
import uuid
import warnings
from datetime import datetime, timedelta

import anndata
import numpy as np
import pandas as pd
import uvicorn
from fastapi import Depends, FastAPI, File, Form, HTTPException, UploadFile
from google.cloud import bigquery

from casp.services import settings
from casp.services.api import async_client, data_controller, exceptions, matching_engine_client, schemas
from casp.services.api.auth import authenticate_user
from casp.services.db import init_db, models, ops

app = FastAPI()
db_session = init_db()


def __log(s):
    dt_string = datetime.now().strftime("%Y/%m/%d %H:%M:%S")
    print(f"{dt_string} - {s}")


async def __get_embeddings(file: bytes, model_name: str) -> t.Tuple[t.List[str], np.array]:
    r_text = await async_client.CASAPIAsyncClient.call_model_service(file_to_embed=file, model_name=model_name)
    # TODO: json isn't very efficient for this
    df = pd.read_json(r_text)
    query_ids = df.obs_ids.tolist()
    embeddings = df.iloc[:, 1:].to_numpy(dtype=float)

    return query_ids, embeddings


def __get_appx_knn_matches(embeddings, model_name: str):
    cas_model = ops.get_model_by(model_name=model_name)
    cas_matching_engine_index = cas_model.cas_matching_engine

    index_endpoint = matching_engine_client.MatchingEngineIndexEndpointGRPCOptionsExposed(
        index_endpoint_name=cas_matching_engine_index.endpoint_id
    )

    response = index_endpoint.match(
        deployed_index_id=cas_matching_engine_index.deployed_index_id,
        queries=embeddings,
        num_neighbors=cas_matching_engine_index.num_neighbors,
    )
    return response


def __get_cell_type_distribution(
    query_ids: t.List, knn_response: t.List, model_name: str, include_dev_metadata: bool
) -> t.List[t.Dict[str, t.Any]]:
    bq_client = bigquery.Client()
    cas_model = ops.get_model_by(model_name=model_name)
    # create temporary table
    my_uuid = str(uuid.uuid4())[:8]
    temp_table_fqn = f"{cas_model.bq_temp_table_dataset}.api_request_{my_uuid}"

    __log(f"Creating Temporary Table {temp_table_fqn}")
    schema = [
        bigquery.SchemaField("query_id", "STRING", mode="REQUIRED"),
        bigquery.SchemaField("match_cas_cell_index", "INT64", mode="REQUIRED"),
        bigquery.SchemaField("match_score", "FLOAT", mode="REQUIRED"),
    ]

    table = bigquery.Table(temp_table_fqn, schema=schema)
    table.expires = datetime.now() + timedelta(minutes=15)
    bq_client.create_table(table)  # Make an API request.

    __log(f"Loading Data Into {temp_table_fqn}")
    __log(f"kNN response length {len(knn_response)}")

    rows_to_insert = []
    for i in range(0, len(knn_response)):
        query_id = query_ids[i]
        for match in knn_response[i]:
            rows_to_insert.append(tuple([str(query_id), int(match.id), float(match.distance)]))

    df = pd.DataFrame(rows_to_insert, columns=["query_id", "match_cas_cell_index", "match_score"])

    job_config = bigquery.LoadJobConfig(schema=schema)
    job = bq_client.load_table_from_dataframe(df, temp_table_fqn, job_config=job_config)
    job.result()  # Wait for the job to complete.

    __log("Querying Match Cell Metadata")
    if include_dev_metadata:
        results = data_controller.get_match_query_metadata_dev_details(
            cas_model=cas_model, match_temp_table_fqn=temp_table_fqn
        )
    else:
        results = data_controller.get_match_query_metadata(cas_model=cas_model, match_temp_table_fqn=temp_table_fqn)

    return results


async def __annotate(
    file: t.BinaryIO, model_name: str, include_dev_metadata: bool = False
) -> t.List[t.Dict[str, t.Any]]:
    file_bytes = io.BytesIO(file.read())
    adata = anndata.read_h5ad(file_bytes)
    file_bytes.seek(0)
    adata_length = adata.shape[0]

    __log("Calculating Embeddings")
    query_ids, embeddings = await __get_embeddings(file_bytes, model_name=model_name)
    if len(query_ids) != adata_length:
        warnings.warn(
            f"Embedding service returned a wrong number of cells. Expected {adata_length}, got {len(query_ids)}"
        )
    __log("Done")
    __log("Performing kNN lookup")
    knn_response = __get_appx_knn_matches(embeddings, model_name=model_name)
    if len(knn_response) != adata_length:
        warnings.warn(
            f"Matching Engine service returned a wrong number of cells. Expected {adata_length}, "
            f"got {len(knn_response)}"
        )
    __log("Getting Cell Type Distributions")
    d = __get_cell_type_distribution(
        query_ids, knn_response, model_name=model_name, include_dev_metadata=include_dev_metadata
    )
    if len(d) != adata_length:
        warnings.warn(
            f"Cell Type Distribution was returned for a wrong number of cells. Expected {adata_length}, got {len(d)}"
        )
    __log("Finished")
    return d


@app.get("/list-models", response_model=t.List[schemas.CASModel])
async def list_models(
    request_user: models.User = Depends(authenticate_user),
):
    return ops.get_models_for_user(user=request_user)


@app.post("/annotate", response_model=t.Union[t.List[schemas.QueryCellDevDetails], t.List[schemas.QueryCell]])
async def annotate(
    file: UploadFile = File(),
    number_of_cells: int = Form(),
    model_name: str = Form(),
    include_dev_metadata: bool = Form(),
    request_user: models.User = Depends(authenticate_user),
):
    """
    Annotate a single anndata file with Cellarium CAS. Input file should be validated and sanitized according to the
    model schema.

    :param file: Byte object of :class:`anndata.AnnData` file to annotate.
    :param number_of_cells: Number of cells in the input file.
    :param model_name: Model name to use for annotation. See `/list-models` endpoint for available models.
    :param include_dev_metadata: Boolean flag indicating whether to include dev metadata in the response.
    :param request_user: Authorized user object obtained  by token from `Bearer` header.
    :return: JSON response with annotations.
    """
    ops.increment_user_cells_processed(request_user, number_of_cells=number_of_cells)
    try:
        return await __annotate(file.file, model_name=model_name, include_dev_metadata=include_dev_metadata)
    except (exceptions.AnnotationServiceException, exceptions.ServiceAPIException) as e:
        raise HTTPException(status_code=500, detail=str(e))


@app.get("/validate-token")
async def validate_token(_: models.User = Depends(authenticate_user)):
    """
    Validate authorization token from `Bearer` header
    :return: Success message if token is valid, otherwise
        return 401 Unauthorized status code if token is invalid or missing
    """
    return {"response": "Success"}


@app.get("/application-info", response_model=schemas.ApplicationInfo)
async def application_info(_: models.User = Depends(authenticate_user)):
    """
    Get Cellarium CAS application info such as version, default feature schema name.
    """
    return data_controller.get_application_info()


@app.get("/feature-schemas", response_model=t.List[schemas.FeatureSchemaInfo])
async def get_feature_schemas(_: models.User = Depends(authenticate_user)):
    """
    Get list of all Cellarium CAS feature schemas

    :return: List of feature schema info objects
    """
    return data_controller.get_feature_schemas()


@app.get("/feature-schema/{schema_name}", response_model=t.List[str])
async def get_feature_schema_by(schema_name: str):
    """
    Get a specific feature schema by its unique name
    :param schema_name: unique feature schema name
    :return: List of features in a correct order
    """
    return data_controller.get_feature_schema_by(schema_name=schema_name)


if __name__ == "__main__":
    uvicorn.run("server:app", host=settings.SERVER_HOST, port=settings.SERVER_PORT, workers=multiprocessing.cpu_count())

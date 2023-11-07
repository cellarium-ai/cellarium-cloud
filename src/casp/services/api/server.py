import multiprocessing
import typing as t
import uuid
from datetime import datetime, timedelta

import pandas as pd
import uvicorn
from fastapi import Depends, FastAPI, File, Form, UploadFile
from google.cloud import bigquery

from casp.services import settings
from casp.services.api import async_client, data_controller, matching_engine_client, schemas
from casp.services.api.auth import authenticate_user
from casp.services.db import init_db, models, ops

if t.TYPE_CHECKING:
    import numpy

app = FastAPI()
db_session = init_db()


def __log(s):
    dt_string = datetime.now().strftime("%Y/%m/%d %H:%M:%S")
    print(f"{dt_string} - {s}")


async def __get_embeddings(myfile, model_name: str) -> t.Tuple[t.List, "numpy.array"]:
    r_text = await async_client.CASAPIAsyncClient.call_model_service(file_to_embed=myfile.read(), model_name=model_name)
    # TODO: json isn't very efficient for this
    df = pd.read_json(r_text)
    query_ids = df.db_ids.tolist()
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


def __get_cell_type_distribution(query_ids, knn_response, model_name: str):
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
    query = f"""
                SELECT t.query_id,
                       ci.cell_type,
                       MIN(t.match_score) min_distance,
                       MAX(t.match_score) max_distance,
                       APPROX_QUANTILES(t.match_score, 100)[SAFE_ORDINAL(25)] as p25_distance,
                       APPROX_QUANTILES(t.match_score, 100)[SAFE_ORDINAL(50)] as median_distance,
                       APPROX_QUANTILES(t.match_score, 100)[SAFE_ORDINAL(75)] as p75_distance,
                       COUNT(*) cell_count
                FROM `{temp_table_fqn}` t
                JOIN `{cas_model.bq_cell_info_table_fqn}` ci ON t.match_cas_cell_index = ci.cas_cell_index
                GROUP BY 1,2
                ORDER BY 1, 8 DESC
                """

    query_job = bq_client.query(query)

    # TODO: use a StreamingResponse and yield/generator pattern here to avoid building entire response in memory, also compress
    # Stream results back in JSONL format (https://jsonlines.org/)
    results = []

    last_query_id = None
    data = {}
    for row in query_job:
        if last_query_id is None or last_query_id != row["query_id"]:
            # emit data and reset state if this isn't the first time through
            if last_query_id is not None:
                results.append(data)

            data = {}
            data["query_cell_id"] = row["query_id"]
            data["matches"] = []
            last_query_id = data["query_cell_id"]

        x = {}
        x["cell_type"] = row["cell_type"]
        x["cell_count"] = row["cell_count"]
        x["min_distance"] = row["min_distance"]
        x["p25_distance"] = row["p25_distance"]
        x["median_distance"] = row["median_distance"]
        x["p75_distance"] = row["p75_distance"]
        x["max_distance"] = row["max_distance"]
        data["matches"].append(x)

    # handle final output
    results.append(data)

    # and return
    return results


async def __annotate(file, model_name: str):
    __log("Calculating Embeddings")
    query_ids, embeddings = await __get_embeddings(file, model_name=model_name)
    __log("Done")

    __log("Performing kNN lookup")
    knn_response = __get_appx_knn_matches(embeddings, model_name=model_name)

    __log("Getting Cell Type Distributions")
    d = __get_cell_type_distribution(query_ids, knn_response, model_name=model_name)

    __log("Finished")
    return d


@app.get("/list-models", response_model=t.List[schemas.CASModel])
async def list_models(
    request_user: models.User = Depends(authenticate_user),
):
    return ops.get_models_for_user(user=request_user)


@app.post("/annotate", response_model=t.List[schemas.QueryCell])
async def annotate(
    myfile: UploadFile = File(),
    number_of_cells: int = Form(),
    model_name: str = Form(),
    request_user: models.User = Depends(authenticate_user),
):
    ops.increment_user_cells_processed(request_user, number_of_cells=number_of_cells)
    return await __annotate(myfile.file, model_name=model_name)


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
    :param _:
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

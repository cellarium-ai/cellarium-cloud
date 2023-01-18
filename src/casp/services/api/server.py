import multiprocessing
import typing as t
import uuid
from datetime import datetime, timedelta

import pandas as pd
import uvicorn
from fastapi import FastAPI, UploadFile
from google.cloud import aiplatform, bigquery

from casp.services import settings
from casp.services.api import async_client, schemas

if t.TYPE_CHECKING:
    import numpy

app = FastAPI()


def __log(s):
    dt_string = datetime.now().strftime("%Y/%m/%d %H:%M:%S")
    print(f"{dt_string} - {s}")


async def __get_embeddings(myfile) -> t.Tuple[t.List, "numpy.array"]:
    r_text = await async_client.CASAPIAsyncClient.call_model_service(file_to_embed=myfile.read())
    # TODO: json isn't very efficient for this
    df = pd.read_json(r_text)
    query_ids = df.db_ids.tolist()
    embeddings = df.iloc[:, 1:].to_numpy(dtype=float)

    return query_ids, embeddings


def __get_appx_knn_matches(embeddings):
    index_endpoint = aiplatform.MatchingEngineIndexEndpoint(index_endpoint_name=settings.KNN_SEARCH_ENDPOINT_ID)

    response = index_endpoint.match(
        deployed_index_id=settings.KNN_SEARCH_DEPLOYED_INDEX_ID,
        queries=embeddings,
        num_neighbors=settings.KNN_SEARCH_NUM_MATCHES,
    )
    return response


def __get_cell_type_distribution(query_ids, knn_response):
    bq_client = bigquery.Client()

    # create temporary table
    my_uuid = str(uuid.uuid4())[:8]
    temp_table_fqn = f"{settings.BQ_TEMP_TABLE_DATASET}.api_request_{my_uuid}"

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
                JOIN `{settings.BQ_CELL_INFO_TABLE_FQN}` ci ON t.match_cas_cell_index = ci.cas_cell_index
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


async def __annotate(file):
    __log("Calculating Embeddings")
    query_ids, embeddings = await __get_embeddings(file)
    __log("Done")

    __log("Performing kNN lookup")
    knn_response = __get_appx_knn_matches(embeddings)

    __log("Getting Cell Type Distributions")
    d = __get_cell_type_distribution(query_ids, knn_response)

    __log("Finished")
    return d


@app.get("/")
async def root() -> str:
    return "Hello world"


@app.post("/annotate", response_model=t.List[schemas.QueryCell])
async def annotate(myfile: UploadFile):
    return await __annotate(myfile.file)


if __name__ == "__main__":
    uvicorn.run(
        "server:app", host=settings.SERVER_HOST, port=settings.SERVER_PORT, workers=multiprocessing.cpu_count() * 2 + 1
    )

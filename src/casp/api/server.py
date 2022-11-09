import typing as t
import time

import uuid
import pandas as pd
import requests
import uvicorn
from fastapi import FastAPI, Form, UploadFile
from google.cloud import aiplatform, bigquery
from pydantic import BaseSettings
from datetime import datetime, timedelta


# TODO --refactor packages, move into commons settings class, leverage .env file
class Settings(BaseSettings):
    server_host: str = "0.0.0.0"
    server_port: int = 8000
    model_server_url: str = "https://casp-pca-serving-vi7nxpvk7a-uk.a.run.app/predict"
    knn_search_endpoint_id: str = "projects/350868384795/locations/us-central1/indexEndpoints/2348891088464379904"
    knn_search_deployed_index_id: str = "deployed_4m_casp_index_v1"
    knn_search_num_matches: int = 100
    bq_cell_info_table_fqn: str = "dsp-cell-annotation-service.cas_4m_dataset.cas_cell_info"
    bq_temp_table_dataset: str = "dsp-cell-annotation-service.cas_4m_dataset" 
    items_per_user: int = 50


settings = Settings()
app = FastAPI()
    
def __log(s):
    dt_string = datetime.now().strftime("%Y/%m/%d %H:%M:%S")
    print(f"{dt_string} - {s}")
    
def __get_embeddings(myfile):
    r = requests.post(settings.model_server_url, files={"file": myfile})

    # this is actually a pandas dataframe in json
    # TODO: json isn't very efficient for this
    df = pd.read_json(r.text)

    query_ids = df.db_ids.tolist()
    embeddings = df.iloc[:, 1:].to_numpy(dtype=float)

    return (query_ids, embeddings)


def __get_appx_knn_matches(embeddings):
    index_endpoint = aiplatform.MatchingEngineIndexEndpoint(index_endpoint_name=settings.knn_search_endpoint_id)

    response = index_endpoint.match(
        deployed_index_id=settings.knn_search_deployed_index_id,
        queries=embeddings,
        num_neighbors=settings.knn_search_num_matches,
    )
    return response


def __get_cell_type_distribution(query_ids, knn_response):
    bq_client = bigquery.Client()

    # create temporary table
    my_uuid = str(uuid.uuid4())[:8]
    temp_table_fqn = f"{settings.bq_temp_table_dataset}.api_request_{my_uuid}"

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

    df = pd.DataFrame(rows_to_insert, columns =['query_id', 'match_cas_cell_index', 'match_score'])

    job_config = bigquery.LoadJobConfig(schema=schema)
    job = bq_client.load_table_from_dataframe(df, temp_table_fqn, job_config=job_config)  
    job.result()  # Wait for the job to complete.

    __log(f"Querying Match Cell Metadata")
    query =    f"""
                SELECT t.query_id, ci.cell_type, avg(t.match_score), count(*) cell_count
                FROM `{temp_table_fqn}` t
                JOIN `{settings.bq_cell_info_table_fqn}` ci ON t.match_cas_cell_index = ci.cas_cell_index
                GROUP BY 1,2
                ORDER BY 1, 3 DESC
                """


    query_job = bq_client.query(query)

    results = {}
    for row in query_job:
        query_id = row["query_id"]
        if query_id not in results:
            results[query_id] = {}

        results[query_id][row["cell_type"]] = row["cell_count"]

    return results


def __annotate(file):
    __log("Calculating Embeddings")
    (query_ids, embeddings) = __get_embeddings(file)

    __log("Performing kNN lookup")
    knn_response = __get_appx_knn_matches(embeddings)

    __log("Getting Cell Type Distributions")
    d = __get_cell_type_distribution(query_ids, knn_response)

    __log("Finished")
    return d


@app.post("/annotate")
async def annotate(myfile: UploadFile, json: str = Form()) -> t.Dict:
    return __annotate(myfile.file)


if __name__ == "__main__":
    uvicorn.run("server:app", host=settings.server_host, port=settings.server_port)

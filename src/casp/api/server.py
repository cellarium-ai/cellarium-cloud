import typing as t

import pandas as pd
import requests
import uvicorn
from fastapi import FastAPI, Form, UploadFile
from google.cloud import aiplatform, bigquery
from pydantic import BaseSettings


# TODO --refactor packages, move into commons settings class, leverage .env file
class Settings(BaseSettings):
    server_host: str = "0.0.0.0"
    server_port: int = 8000
    model_server_url: str= "https://casp-pca-serving-vi7nxpvk7a-uk.a.run.app/predict"
    knn_search_endpoint_id: str = "projects/350868384795/locations/us-central1/indexEndpoints/2348891088464379904"
    knn_search_deployed_index_id: str = "deployed_4m_casp_index_v1"
    knn_search_num_matches: int = 100
    bq_cell_info_table_fqn: str = "dsp-cell-annotation-service.cas_4m_dataset.cas_cell_info"
    items_per_user: int = 50

settings = Settings()
app = FastAPI()

def __get_embeddings(myfile):
    r = requests.post(settings.model_server_url, files={'file': myfile})

    # this is actually a pandas dataframe in json
    # TODO: json isn't very efficient for this    
    df = pd.read_json(r.text)

    query_ids = df.db_ids.tolist()
    embeddings = df.iloc[:, 1:].to_numpy(dtype=float)

    return (query_ids, embeddings)

def __get_appx_knn_matches(embeddings):
    index_endpoint = aiplatform.MatchingEngineIndexEndpoint(index_endpoint_name=settings.knn_search_endpoint_id)

    response = index_endpoint.match(deployed_index_id=settings.knn_search_deployed_index_id, 
                                    queries=embeddings, 
                                    num_neighbors=settings.knn_search_num_matches)
    return response

# TODO: instead upload query_id, match_cas_cell_id into a temporary table
# and then join to that rather than batch these queries.
def __get_cell_type_distribution(query_ids, knn_response):
    bq_client = bigquery.Client()  
    query_pieces = []

    for i in range(0,len(knn_response)):
        query_id = query_ids[i]
        result=knn_response[i]
    
        match_ids = [match.id for match in result]

        query_pieces.append(f"""
            SELECT '{query_id}' as query_id, cell_type, count(*) cell_count
            FROM {settings.bq_cell_info_table_fqn}
            WHERE cas_cell_index IN ({",".join(match_ids)})
            GROUP BY 1,2
         """
         )
    

    results = {}
    batch_size = 100
    for i in range(0, len(query_pieces), batch_size):
        print(f"Retrieving chunk {int(i/batch_size)} from BigQuery")
        chunk = query_pieces[i:i + batch_size]
        final_query = " UNION ALL ".join(chunk) + " ORDER BY 1, 3 DESC"
        query_job = bq_client.query(final_query)

        for row in query_job:
            query_id = row["query_id"]
            if query_id not in results:
                results[query_id] = {}
        
            results[query_id][row["cell_type"]] = row["cell_count"]

    return results

def __annotate(file):
    print("TODO: Validation and Preprocessing")
    #TODO get allow-list of ensembl ids and subset

    print("Calculating Embeddings")
    (query_ids, embeddings) = __get_embeddings(file)

    print("Performing kNN lookup")
    knn_response = __get_appx_knn_matches(embeddings)

    print("Getting Cell Type Distributions")
    d = __get_cell_type_distribution(query_ids, knn_response)
    return d


@app.post("/annotate")
async def annotate(myfile: UploadFile, json: str = Form()) -> t.Dict:
    return __annotate(myfile.file)

if __name__ == "__main__":
    uvicorn.run("server:app", host=settings.server_host, port=settings.server_port)

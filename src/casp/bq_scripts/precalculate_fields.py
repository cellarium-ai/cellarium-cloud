import typing as t

from google.cloud import bigquery
from dask_bigquery import read_gbq
from pandas_gbq import to_gbq
import numpy as np
import json
import time

from casp.bq_scripts import constants

SQL_PRECALCULATE_TOTAL_MRNA_UMIS_FORMAT = f"""
    UPDATE `{{project}}.{{dataset}}.cas_cell_info` ci
    SET ci.total_mrna_umis = m.raw_counts_total
    FROM (
        SELECT cas_cell_index, SUM(raw_counts) AS raw_counts_total
        FROM `{{project}}.{{dataset}}.cas_raw_count_matrix` m
        JOIN `{{project}}.{{dataset}}.cas_feature_info` fi ON (fi.cas_feature_index = m.cas_feature_index)
        WHERE fi.feature_biotype = '{constants.MRNA_FEATURE_BIOTYPE_NAME}'
        GROUP BY cas_cell_index
    ) AS m
    WHERE ci.cas_cell_index = m.cas_cell_index;
"""

SQL_FIELD_MAPPING = {"total_mrna_umis": SQL_PRECALCULATE_TOTAL_MRNA_UMIS_FORMAT}


def precalculate_batch(client: bigquery.Client, dataset: str, project: str):
    # since I cannot figure out how to do this with one query, I will do two and use pandas locally and then upload
    ingest_query_job = client.query(
        f"SELECT cas_ingest_id, uns_metadata FROM `{project}.{dataset}.cas_ingest_info`"
    )
    print("Querying ingest info")
    df_ingest = ingest_query_job.to_dataframe()

    # annotate batch_key per ingest
    df_ingest['batch_key'] = df_ingest['uns_metadata'].apply(
        lambda x: json.loads(x).get('batch_condition', ['donor_id'])[0],
    )

    destination_table = f'{project}.{dataset}.tmp__cell_info'
    try:
        # read columns of obs to a temporary table and pull that table in full from bigquery
        job_config = bigquery.QueryJobConfig(destination=destination_table)
        obs_query_job = client.query(
            f"""SELECT cas_cell_index, cas_ingest_id, donor_id, obs_metadata_extra 
            FROM `{project}.{dataset}.cas_cell_info` ci
            WHERE ci.is_primary_data""",
            job_config=job_config,
        )
        print("Querying cell metadata")
        obs_query_job.result()
        df_obs = read_gbq(project_id=project, dataset_id=dataset, table_id='cas_cell_info')
    except:
        pass
    client.delete_table(destination_table)

    # merge the two dataframes
    df = df_obs.merge(df_ingest, on='cas_ingest_id', how='left')

    # compute the batch labels for each cell
    df['batch'] = df.apply(
        lambda d: (
            d['cas_ingest_id'] + '::' 
            + json.loads(d['obs_metadata_extra']).get(d['batch_key'], d['donor_id'])
        ), 
        axis=1, 
        meta=str,  # dask
    )

    def next_int():
        i = 0
        while True:
            yield i
            i += 1

    get_next_int = next_int()

    # upload a table with cas_cell_index and batch
    def upload_chunk(chunk):
        time.sleep(next(get_next_int) * 2)  # avoid hitting bigquery api rate limits
        to_gbq(chunk, destination_table=f'{dataset}.tmp__cell_info_batch', project_id=project, if_exists='append')

    df = df.repartition(npartitions=20)  # ensure we don't hit a bigquery api rate limit with too many partitions
    df[['cas_cell_index', 'batch']].map_partitions(upload_chunk, meta=df._meta).compute()

def precalculate_fields(dataset: str, fields: t.List[str], project: str):
    if not set(fields).issubset(set(SQL_FIELD_MAPPING.keys())):
        raise ValueError("Fields should be one of the `SQL_FIELD_MAPPING` keys")

    client = bigquery.Client(project=project)

    for field in fields:
        print(f"Executing calculation for {field}")

        match field:
            case "batch":
                precalculate_batch(client, dataset, project)
                
            case _:
                sql_format = SQL_FIELD_MAPPING[field]
                query_job = client.query(query=sql_format.format(project=project, dataset=dataset))
                _ = query_job.result()  # Wait till command runs

        print("Done.")

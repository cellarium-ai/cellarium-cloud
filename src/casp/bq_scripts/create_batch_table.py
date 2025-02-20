from google.cloud import bigquery
from dask_bigquery import read_gbq
from pandas_gbq import to_gbq
import numpy as np
import json
import time

project = 'dsp-cell-annotation-service'
dataset = 'cas_50m_dataset'
destination_table = f'{dataset}.tmp__cell_info_batch'

client = bigquery.Client(project=project)

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

tmp_destination_table_name = 'tmp__cell_info'
tmp_destination_table = f'{project}.{dataset}.{tmp_destination_table_name}'
# client.delete_table(tmp_destination_table, not_found_ok=True)
# # read columns of obs to a temporary table and pull that table in full from bigquery
# job_config = bigquery.QueryJobConfig(destination=tmp_destination_table)
# obs_query_job = client.query(
#     f"""SELECT cas_cell_index, cas_ingest_id, donor_id, obs_metadata_extra 
#     FROM `{project}.{dataset}.cas_cell_info` ci
#     WHERE ci.is_primary_data""",
#     job_config=job_config,
# )
# print("Querying cell metadata")
# obs_query_job.result()
df_obs = read_gbq(project_id=project, dataset_id=dataset, table_id=tmp_destination_table_name)


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

client.delete_table(destination_table, not_found_ok=True)

# upload a table with cas_cell_index and batch
def upload_chunk(chunk):
    time.sleep(next(get_next_int) * 1)  # avoid hitting bigquery api rate limits
    to_gbq(chunk, destination_table=destination_table, project_id=project, if_exists='append')

print(f"Uploading to {destination_table}")
df = df.repartition(npartitions=1000)  # ensure we don't hit a bigquery api rate limit with too many partitions
df[['cas_cell_index', 'batch']].map_partitions(upload_chunk, meta=df._meta).compute()

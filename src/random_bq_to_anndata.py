from google.cloud import bigquery
import argparse
import random
import anndata as ad


def get_random_ids(project, dataset, num_cells):
    client = bigquery.Client()
    query = client.query(f"SELECT MAX(cas_cell_index) AS max_table_number FROM `{project}.{dataset}.cas_cell_info`")
    max_cell_id = int([row.max_table_number for row in list(query.result())][0])

    print(f"Getting {num_cells} random IDs between 0 and {max_cell_id}...")
    return random.choices(range(0, max_cell_id, 1), k = num_cells)


def get_cell_data(project, dataset, num_cells):
    random_ids = get_random_ids(project, dataset, num_cells)

    in_clause = f"WHERE cas_cell_index IN ({','.join(map(str, random_ids))})"
    print(in_clause)
    client = bigquery.Client()
    # query = client.query(f"SELECT ? FROM `{project}.{dataset}.cas_cell_info` " + in_clause)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(allow_abbrev=False, description='Query CASP tables for random cells')
    parser.add_argument('--project', type=str, help='BigQuery Project', required=True)
    parser.add_argument('--dataset', type=str, help='BigQuery Dataset', required=True)
    parser.add_argument('--num_cells', type=int, help='Number of cells to return', required=True)

    args = parser.parse_args()
    get_cell_data(args.project, args.dataset, args.num_cells)

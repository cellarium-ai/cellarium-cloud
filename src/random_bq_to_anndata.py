from google.cloud import bigquery
import argparse
import random
import anndata as ad


def get_random_ids(project, dataset, client, num_cells):
    query = client.query(f"SELECT MAX(cas_cell_index) AS max_table_number FROM `{project}.{dataset}.cas_cell_info`")
    max_cell_id = int([row.max_table_number for row in list(query.result())][0])
    print(f"Getting {num_cells} random IDs between 0 and {max_cell_id}...")
    cell_ids = list(range(0, max_cell_id + 1))
    random.shuffle(cell_ids)
    del cell_ids[num_cells:]
    return cell_ids


def get_cell_data(project, dataset, client, num_cells):
    random_ids = get_random_ids(project, dataset, client, num_cells)
    in_clause = f" matrix.cas_cell_index IN ({','.join(map(str, random_ids))})"

    # at some point, we will probably want create temp table of cell_ids and then JOIN on it
    # instead of an IN clause
    sql = f"SELECT original_cell_id, cell_type, original_gene_id, feature_name AS gene_feature, raw_counts AS count FROM `{project}.{dataset}.cas_cell_info` AS cell, `{project}.{dataset}.cas_gene_info` AS gene, `{project}.{dataset}.cas_raw_count_matrix` AS matrix WHERE matrix.cas_cell_index = cell.cas_cell_index AND matrix.cas_gene_index = gene.cas_gene_index AND" + in_clause
    print(f"Getting {num_cells} random cells' data from {project}.{dataset}...")
    query = client.query(sql)
    return query.result()


def random_bq_to_anndata(project, dataset, num_cells):
    client = bigquery.Client(project=project)
    cell_data = get_cell_data(project, dataset, client, num_cells)
    adata = ad.AnnData()
    for row in list(cell_data):
        ...


if __name__ == '__main__':
    parser = argparse.ArgumentParser(allow_abbrev=False, description='Query CASP tables for random cells')
    parser.add_argument('--project', type=str, help='BigQuery Project', required=True)
    parser.add_argument('--dataset', type=str, help='BigQuery Dataset', required=True)
    parser.add_argument('--num_cells', type=int, help='Number of cells to return', required=True)

    args = parser.parse_args()
    random_bq_to_anndata(args.project, args.dataset, args.num_cells)

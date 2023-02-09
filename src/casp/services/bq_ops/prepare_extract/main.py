import argparse

from casp.bq_scripts import prepare_extract
from casp.services import utils

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        allow_abbrev=False, description="Prepare CASP tables ML Training/Inference Extract"
    )
    parser.add_argument("--dataset", type=str, help="BigQuery Dataset", required=True)
    parser.add_argument("--extract_table_prefix", type=str, help="Prefix for extract tables", required=True)
    parser.add_argument(
        "--min_observed_cells", type=int, help="minimum observed cells per gene", default=3, required=False
    )
    parser.add_argument(
        "--extract_bin_size", type=int, help="desired cells per extract bin", default=10000, required=False
    )
    parser.add_argument(
        "--fq_allowed_original_feature_ids",
        type=str,
        help="fully qualified reference to table of allowed feature names",
        default="dsp-cell-annotation-service.cas_reference_data.refdata-gex-GRCh38-2020-A",
        required=False,
    )
    credentials, project_id = utils.get_google_service_credentials()
    args = parser.parse_args()

    prepare_extract(
        project=project_id,
        dataset=args.dataset,
        extract_table_prefix=args.extract_table_prefix,
        min_observed_cells=args.min_observed_cells,
        fq_allowed_original_feature_ids=args.fq_allowed_original_feature_ids,
        extract_bin_size=args.extract_bin_size,
    )

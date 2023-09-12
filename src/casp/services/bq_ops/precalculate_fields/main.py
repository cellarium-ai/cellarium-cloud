import argparse

from casp.bq_scripts import precalculate_fields


def main(dataset: str, fields_str: str):
    fields_list = fields_str.split(",")
    precalculate_fields(dataset=dataset, fields=fields_list)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(allow_abbrev=False)

    parser.add_argument("--dataset", type=str, help="BigQuery Dataset", required=True)
    parser.add_argument(
        "--fields_str",
        type=str,
        help="Fields to precalculate. Should be one of `SQL_FIELD_MAPPING` keys. Comma separated",
        required=True,
    )
    args = parser.parse_args()

    main(dataset=args.dataset, fields_str=args.fields_str)

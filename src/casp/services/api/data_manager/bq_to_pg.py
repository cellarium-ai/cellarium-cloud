import argparse
import concurrent.futures
import json
import logging
import time
import traceback
import typing as t
from uuid import uuid4

from google.cloud import bigquery_storage
from google.cloud.bigquery_storage import ReadStream, types

from casp.data_manager.base_data_manager import BaseDataManager
from casp.services import db, settings
from casp.services.db import models
from casp.services.utils import get_google_service_credentials

BIGQUERY_TABLE = "bigquery_table"
POSTGRES_TABLE = "postgres_table"
COLUMNS_STRUCTURE = "cols"
CHUNK_SIZE = "chunk_size"


logger = logging.getLogger(__name__)
logging.basicConfig(level=settings.LOG_LEVEL.upper())


class BQToPG(BaseDataManager):
    """
    Copy data from BigQuery to PostgreSQL in chunks.

    This is very loosely based off of the Medium article:
    https://ysrazsingh.medium.com/efficient-data-transfer-migrating-from-bigquery-to-postgresql-90a3f70127bd
    """

    def _read_json_data(self, json_file_path):
        """Read JSON migration data from a file."""
        with open(json_file_path) as f:
            return json.load(f)

    def _write_chunk(self, table: db.Base, data: t.List[t.Dict[str, any]]) -> None:
        with self.system_data_db_session_maker() as session:
            session.bulk_insert_mappings(table, data)
            session.commit()

    def _process_stream(
        self,
        reader: ReadStream,
        stream_num: int,
        stream_name: str,
        table: db.Base,
        columns,
        ingest_id: str,
        chunk_size: int,
    ) -> None:
        logger.info(f"Processing stream {stream_num}:{stream_name} with chunk size of {chunk_size}")
        with concurrent.futures.ThreadPoolExecutor(1) as executor:
            # This really will only ever have one future at any given time
            futures = []
            data = []
            try:
                for i, row in enumerate(reader.rows()):
                    if i % chunk_size == 0 and i != 0:
                        # See if any previous writes failed
                        for future in futures:
                            future.result()
                        futures.append(executor.submit(self._write_chunk, table=table, data=data))
                        logger.info(f"Stream {stream_num}:{stream_name}: Wrote {i} records")
                        data = []
                    data_dict = {"cas_ingest_id": ingest_id}
                    data_dict["cas_cell_index"] = i
                    for col in columns:
                        val = row[col["source_column"]].as_py()
                        col_value = col.get("destination_column", col["source_column"])
                        data_dict[col_value] = str(val)
                    data.append(data_dict)

                if len(data) > 0:
                    # See if any previous writes failed
                    for future in futures:
                        future.result()
                    futures.append(executor.submit(self._write_chunk, table=table, data=data))
                    logger.info(f"Stream {stream_num}:{stream_name}: Wrote {i} records")
                concurrent.futures.wait(futures, return_when=concurrent.futures.ALL_COMPLETED)

                # See if any previous writes failed
                for future in futures:
                    future.result()
            except Exception as e:
                traceback.print_exc()
                logger.error(f"Error processing stream {stream_name}")
                raise e

    def _get_table(self, table_name: str) -> db.Base:
        tables = db.Base.__subclasses__()
        for table in tables:
            if table.__tablename__ == table_name:
                logger.info(f"Writing to table {table_name}")
                return table
        raise Exception(f"Table {table_name} is not recognized")

    def ingest(self, blueprint: str = settings.COPY_BLUEPRINT, concurrency: int = settings.COPY_NUM_STREAMS):
        blueprint = blueprint or settings.COPY_BLUEPRINT
        concurrency = concurrency or settings.COPY_NUM_STREAMS
        logger.info(f"Starting ingest blueprint {blueprint} with max conurrency {concurrency}...")
        start = time.time()
        json_file_path = f"{settings.APP_ROOT}/casp/services/api/data_manager/bq_to_pg_blueprints/{blueprint}.json"

        project_id = get_google_service_credentials()[1]

        json_data = self._read_json_data(json_file_path)

        for data in json_data:
            source_table = data[BIGQUERY_TABLE]
            destination_table_name = data[POSTGRES_TABLE]
            destination_table = self._get_table(destination_table_name)
            columns = data[COLUMNS_STRUCTURE]
            chunk_size = data[CHUNK_SIZE]
            table = f"projects/{project_id}/datasets/{settings.DEFAULT_MODEL_BQ_DATASET_NAME}/tables/{source_table}"

            read_options = types.ReadSession.TableReadOptions(selected_fields=[col["source_column"] for col in columns])

            parent = "projects/{}".format(project_id)

            requested_session = types.ReadSession(
                table=table,
                data_format=types.DataFormat.ARROW,
                read_options=read_options,
            )
            client = bigquery_storage.BigQueryReadClient()
            read_session = client.create_read_session(
                parent=parent,
                read_session=requested_session,
                max_stream_count=concurrency,
            )

            ingest_id = None
            with self.system_data_db_session_maker() as session:
                ingest_info = models.CellIngestInfo(
                    cas_ingest_id=uuid4(),
                    dataset_id=f"{settings.DEFAULT_MODEL_BQ_DATASET_NAME}.{source_table}",
                    uns_metadata={},
                )

                session.add(ingest_info)
                ingest_id = ingest_info.cas_ingest_id
                logger.info(f"Created ingest id {ingest_id}")
                session.commit()

            with concurrent.futures.ThreadPoolExecutor(concurrency) as executor:
                futures = [
                    executor.submit(
                        self._process_stream,
                        reader=client.read_rows(stream.name),
                        stream_num=stream_num,
                        stream_name=stream.name,
                        table=destination_table,
                        columns=columns,
                        ingest_id=ingest_id,
                        chunk_size=chunk_size,
                    )
                    for stream_num, stream in enumerate(read_session.streams)
                ]
                concurrent.futures.wait(futures, return_when=concurrent.futures.ALL_COMPLETED)
                try:
                    # See if any streams failed
                    for future in futures:
                        future.result()
                    logger.info(f"Data copy completed successfully for table {destination_table_name}!")
                except Exception as e:
                    logger.info(f"Data copy failed for table {destination_table_name}!")
                    raise e
            logger.info(f"Time taken: {time.time() - start} seconds")


if __name__ == "__main__":
    bq_to_pg = BQToPG()

    parser = argparse.ArgumentParser()
    parser.add_argument("--concurrency", dest="concurrency", type=int, help="The number of streams to use for the copy")
    parser.add_argument(
        "--blueprint", dest="blueprint", type=str, help="The blueprint (e.g. column mappings) to use for the copy"
    )
    args = parser.parse_args()

    bq_to_pg.ingest(concurrency=args.concurrency, blueprint=args.blueprint)

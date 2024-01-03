from casp.services.api_internal.data_manager import EmbeddingModelRegistryDataManager
from casp.services.db import models


class EmbeddingModelRegistryService:
    def __init__(self):
        self.data_manager = EmbeddingModelRegistryDataManager()

    def create_embedding_model(
        self, model_name: str, model_file_path: str, embedding_dimension: int, bq_dataset_name: str, schema_name: str
    ) -> models.CASModel:
        return self.data_manager.create_embedding_model(
            model_name=model_name,
            model_file_path=model_file_path,
            embedding_dimension=embedding_dimension,
            bq_dataset_name=bq_dataset_name,
            schema_name=schema_name,
        )

    def create_index(
        self,
        model_name: str,
        index_name: str,
        deployed_index_id: str,
        endpoint_id: str,
        embedding_dimension: int,
        num_neighbors: int,
    ) -> models.CASMatchingEngineIndex:
        return self.data_manager.create_space_vector_search_index(
            model_name=model_name,
            index_name=index_name,
            deployed_index_id=deployed_index_id,
            endpoint_id=endpoint_id,
            embedding_dimension=embedding_dimension,
            num_neighbors=num_neighbors,
        )

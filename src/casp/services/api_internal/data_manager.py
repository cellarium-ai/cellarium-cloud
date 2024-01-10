from sqlalchemy.exc import IntegrityError

from casp.data_manager import BaseDataManager
from casp.services.api_internal import exceptions
from casp.services.db import models


class EmbeddingModelRegistryDataManager(BaseDataManager):
    def create_embedding_model(
        self, model_name: str, model_file_path: str, embedding_dimension: int, bq_dataset_name: str, schema_name: str
    ) -> models.CASModel:
        """
        Create a new embedding model in the database

        :param model_name: Model name
        :param model_file_path: Model file path in GCS
        :param embedding_dimension: Model embedding dimension
        :param bq_dataset_name: BigQuery dataset name which was used to make data extract for training the model
        :param schema_name: Schema name which was used to make data extract for training the model
        """
        try:
            new_model = models.CASModel(
                model_name=model_name,
                embedding_dimension=embedding_dimension,
                model_file_path=model_file_path,
                bq_dataset_name=bq_dataset_name,
                schema_name=schema_name,
            )
            self.postgres_db_session.add(new_model)
            self.postgres_db_session.commit()
        except IntegrityError as e:
            raise exceptions.ModelUniqueConstraintError(
                f"{model_name} or {model_file_path} already present in the database"
            ) from e
        return models.CASModel.query.filter_by(model_name=model_name).first()

    def create_space_vector_search_index(
        self,
        model_name: str,
        index_name: str,
        deployed_index_id: str,
        endpoint_id: str,
        embedding_dimension: int,
        num_neighbors: int,
    ) -> models.CASMatchingEngineIndex:
        """
        Create a new space vector search index in the database

        :param model_name: Model name to which the index belongs
        :param index_name: Index name
        :param deployed_index_id: Deployed index ID in GCP Vertex AI
        :param endpoint_id: Endpoint ID in GCP Vertex AI
        :param embedding_dimension: Embedding dimension of the index and the model
        :param num_neighbors: Number of neighbors that will be returned by the index
        """
        model = models.CASModel.query.filter_by(model_name=model_name).first()

        if model is None:
            raise exceptions.ModelNotFoundError(f"Model `{model_name}` not found in the database")

        try:
            new_index = models.CASMatchingEngineIndex(
                index_name=index_name,
                deployed_index_id=deployed_index_id,
                endpoint_id=endpoint_id,
                embedding_dimension=embedding_dimension,
                num_neighbors=num_neighbors,
                model_id=model.id,
            )
            self.postgres_db_session.add(new_index)
            self.postgres_db_session.commit()
        except IntegrityError as e:
            raise exceptions.IndexUniqueConstraintError(f"`{index_name}` already present in the database") from e

        return models.CASMatchingEngineIndex.query.filter_by(index_name=index_name).first()

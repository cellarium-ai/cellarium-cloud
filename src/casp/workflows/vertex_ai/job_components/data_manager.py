from casp.data_manager import BaseDataManager
from casp.services import settings
from casp.services.db import models


class EmbeddingModelRegistryDataManager(BaseDataManager):
    def register_embedding_model(self, model_name: str, model_file_path: str, embedding_dimension: int):
        new_model = models.CASModel(
            model_name=model_name,
            embedding_dimension=embedding_dimension,
            model_file_path=model_file_path,
            admin_use_only=True,
            is_default_model=False,
        )
        self.postgres_db_session.add(new_model)
        self.postgres_db_session.commit()

    def register_space_vector_search_index(
        self, model_name: str, index_name: str, deployed_index_id: str, endpoint_id: str
    ):
        model = models.CASModel.query.filter_by(model_name=model_name).first()
        new_index = models.CASMatchingEngineIndex(
            index_name=index_name,
            deployed_index_id=deployed_index_id,
            endpoint_id=endpoint_id,
            embedding_dimension=model.embedding_dimension,
            num_neighbors=settings.KNN_SEARCH_NUM_MATCHES_DEFAULT,
            model_id=model.id,
        )
        self.postgres_db_session.add(new_index)
        self.postgres_db_session.commit()

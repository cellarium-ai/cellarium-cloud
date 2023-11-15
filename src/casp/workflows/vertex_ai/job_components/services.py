from casp.workflows.vertex_ai.job_components.data_manager import EmbeddingModelRegistryDataManager


class EmbeddingModelRegistryService:
    def __init__(self):
        self.data_manager = EmbeddingModelRegistryDataManager()

    def register_embedding_model(self, model_name: str, model_file_path: str, embedding_dimension: int):
        self.data_manager.register_embedding_model(
            model_name=model_name, model_file_path=model_file_path, embedding_dimension=embedding_dimension
        )

    # def create_space_vector_search_index(self, model_name: str):
    #     project_id, credentials = utils.get_google_service_credentials()
    #     aiplatform.init(project=project_id, credentials=credentials)
    #
    #     index = aiplatform.MatchingEngineIndex.create(
    #         display_name=display_name,
    #         metadata_schema_uri="gs://google-cloud-aiplatform/schema/index/vertex-ai-matching-engine-brute-force-v1.yaml",
    #         metadata={"contentsDeltaUri": gcs_uri}
    #     )
    #
    #     index.wait()

import typing as t

from google.oauth2.service_account import Credentials
import google.cloud.run_v2 as run_v2

from casp.services import settings


def get_google_service_credentials() -> t.Tuple["Credentials", str]:
    credentials = Credentials.from_service_account_info(
        info=settings.GOOGLE_ACCOUNT_CREDENTIALS, scopes=None, default_scopes=None
    )
    return credentials, settings.GOOGLE_ACCOUNT_CREDENTIALS.get("project_id")


def deploy_cloud_run_model(model_file_path):
    # Create a client
    credentials, project_id = get_google_service_credentials()
    client = run_v2.ServicesClient(credentials=credentials)
    # Initialize request argument(s)
    template = run_v2.RevisionTemplate(
        containers=[
            run_v2.Container(
                image=settings.CLOUD_RUN_IMAGE_NAME,
                command=["python"],
                args=["casp/services/model_inference/pca/server.py", "--model_file_path", model_file_path],
                ports=[run_v2.ContainerPort(name="http1", container_port=8000)],
                resources=run_v2.ResourceRequirements(limits={"cpu": "4000m", "memory": "8Gi"})
            )],
        vpc_access=run_v2.VpcAccess(connector=settings.VPC_CONNECTOR_NAME)
    )
    service = run_v2.Service(template=template)
    request = run_v2.CreateServiceRequest(
        parent=f"projects/{project_id}/locations/us-central1",
        service_id="cloud-run-from-admin-test",
        service=service
    )
    operation = client.create_service(request=request)
    print("Waiting for operation to complete...")
    response = operation.result()
    # Handle the response
    print(response)

# def deploy_cloud_run():
#     ServicesAsyncClient.create_service()

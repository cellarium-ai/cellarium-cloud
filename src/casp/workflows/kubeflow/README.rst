KubeFlow Pipelines for Cellarium Cloud
######################################

This module contains the code for the KubeFlow Pipelines for Cellarium Cloud. Currently, the module works with
vertex-ai and AI Platform Pipelines. To use the module, you need to have a GCP project and a GCP service account.

All pipeline tasks are executed in Cellarium Cloud Pipeline docker image.
Task structure:
[input config file GCS path] -> (Task) -> [put output (if any) to GCS path]

Pipelines consist of the tasks that are executed in the order they are defined.
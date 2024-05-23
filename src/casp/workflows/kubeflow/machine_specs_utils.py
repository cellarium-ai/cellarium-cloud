import json
import os
import typing as t

from casp.workflows.kubeflow import constants


def create_empty_machine_specs_json_file_if_not_exists():
    directory = os.path.dirname(constants.MACHINE_SPEC_FILE_PATH)

    os.makedirs(directory, exist_ok=True)

    if not os.path.exists(constants.MACHINE_SPEC_FILE_PATH):
        with open(constants.MACHINE_SPEC_FILE_PATH, "w") as file:
            json.dump({}, file)


def read_machine_specs():
    create_empty_machine_specs_json_file_if_not_exists()
    with open(constants.MACHINE_SPEC_FILE_PATH, "rb") as f:
        return json.loads(f.read())


def update_machine_specs(update_with: t.Dict[str, t.Any]) -> None:
    with open(constants.MACHINE_SPEC_FILE_PATH, "r+") as f:
        current_machine_specs = json.loads(f.read())
        f.seek(0)
        current_machine_specs.update(update_with)
        f.write(json.dumps(current_machine_specs))


def destroy_machine_specs_file() -> None:
    os.remove(constants.MACHINE_SPEC_FILE_PATH)


DUMMY_SPEC = {
    "display_name": "Dummy Display Name",
    "replica_count": 1,
    "machine_type": "dummy-type",
    "accelerator_type": None,
    "accelerator_count": 4,
    "base_image": "dummy-docker.pkg.dev/dummy/dummy/dummy-image:dummy-tag",
}

from kfp.dsl.python_component import PythonComponent

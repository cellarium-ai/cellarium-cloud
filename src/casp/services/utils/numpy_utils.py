import base64
import json

import numpy as np


def numpy_to_base64(array: np.array) -> str:
    """
    Convert numpy array to a JSON string containing shape, dtype, and raw data, then encode the JSON string as base64.

    :param array: Numpy array to convert.
    :return: Base64 string of the JSON-encoded array data, shape, and dtype.
    """
    # Create a JSON object containing the shape, dtype, and raw data of the array
    data = {
        "shape": array.shape,
        "dtype": str(array.dtype),
        "data": array.tobytes().hex(),  # Convert the raw data to a hexadecimal string
    }

    # Convert the dictionary to a JSON string and then to a base64 string
    json_data = json.dumps(data)
    return base64.b64encode(json_data.encode("utf-8")).decode("utf-8")


def base64_to_numpy(base64_str: str) -> np.array:
    """
    Convert base64 string back to numpy array, restoring original shape, dtype, and data.

    :param base64_str: Base64 string to convert.
    :return: Numpy array.
    """
    json_str = base64.b64decode(base64_str.encode("utf-8")).decode("utf-8")

    data = json.loads(json_str)
    array_data = bytes.fromhex(data["data"])  # Convert the hexadecimal string back to bytes

    return np.frombuffer(array_data, dtype=data["dtype"]).reshape(data["shape"])

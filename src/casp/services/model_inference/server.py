import multiprocessing
import typing as t

import pandas as pd
import uvicorn
from fastapi import FastAPI, File

from casp.services.model_inference import utils
from casp.services.db import ops

app = FastAPI()


@app.post("/predict")
async def predict(model_system_name: str, file: bytes = File()) -> t.Dict:
    # Get Model dump file
    model_info = ops.get_model_by(system_name=model_system_name)
    dump_manager = utils.get_dump_manager(model_info.model_file_path)
    model = dump_manager.model
    transform = dump_manager.transform
    # Get data and transform it
    X, db_ids = utils.load_data(file)
    X = transform(X)
    # Embed data
    embeddings = model.transform(X)
    df = pd.DataFrame(embeddings.numpy())
    columns = ["db_ids", *list(df.columns.values)]
    df["db_ids"] = db_ids.numpy().astype(int)
    df = df[columns]
    return df.to_dict()


if __name__ == "__main__":
    # Run model server
    uvicorn.run("server:app", host="0.0.0.0", port=8001, )  # workers=multiprocessing.cpu_count() * 2 + 1)

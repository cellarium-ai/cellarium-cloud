import multiprocessing
import typing as t

import pandas as pd
import uvicorn
from fastapi import FastAPI, File

from casp.ml.services.inference.pca import utils

app = FastAPI()

dump_manager = utils.get_dump_manager()
model = dump_manager.model
transform = dump_manager.transform


@app.get("/")
async def read_root():
    return {"message": "Welcome from the API"}


@app.post("/predict")
async def predict(file: bytes = File()) -> t.Dict:
    X, db_ids = utils.load_data(file)
    X = transform(X)
    embeddings = model.transform(X)
    df = pd.DataFrame(embeddings.numpy())
    columns = ["db_ids", *list(df.columns.values)]
    df["db_ids"] = db_ids.numpy().astype(int)
    df = df[columns]
    return df.to_dict()


if __name__ == "__main__":
    uvicorn.run("server:app", host="0.0.0.0", port=8000, workers=multiprocessing.cpu_count() * 2 + 1)

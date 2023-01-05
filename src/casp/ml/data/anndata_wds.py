import webdataset as wds
from webdataset import filters, gopen
from webdataset.handlers import reraise_exception
from webdataset.tariterators import url_opener, group_by_keys

import anndata
from io import DEFAULT_BUFFER_SIZE
import torch
import torch.utils.data
import torch.nn
import os

#pip install Pillow
#pip install webdataset

os.environ["WDS_VERBOSE_CACHE"] = "1"

def baseline():
    url = "http://storage.googleapis.com/nvdata-publaynet/publaynet-train-{000000..000009}.tar"
    dataset = wds.WebDataset(url).shuffle(1000).decode("rgb").to_tuple("png", "json")
    isinstance(dataset, torch.utils.data.IterableDataset)
    for image, json in dataset:
        break
    print(json)

def use_gs_urls():
    ## Experiment 1: move public http-based tar files to protected google bucket
    # gsutil cp gs://nvdata-publaynet/publaynet-train-00000?.tar gs://dsp-cell-annotation-service/nvdata-publaynet-mirror/
    url = "gs://dsp-cell-annotation-service/nvdata-publaynet-mirror/publaynet-train-{000000..000009}.tar"
    dataset = wds.WebDataset(url).shuffle(1000).decode("rgb").to_tuple("png", "json")
    isinstance(dataset, torch.utils.data.IterableDataset)
    for image, json in dataset:
        break
    print(json)

def use_data_pipeline():
    ## Experiment 2: unroll WebDataset into individual pipeline operations (so we can move to h5ad)
    url = "gs://dsp-cell-annotation-service/nvdata-publaynet-mirror/publaynet-train-{000000..000009}.tar"

    dataset2 = wds.DataPipeline(
        wds.SimpleShardList(url),
        # at this point we have an iterator over all the shards
        wds.shuffle(100),
        wds.split_by_worker,
        # at this point, we have an iterator over the shards assigned to each worker
        wds.tarfile_to_samples(),
        wds.shuffle(1000),
        wds.decode("rgb"),
        wds.to_tuple("png", "json")
    )

    isinstance(dataset2, torch.utils.data.IterableDataset)
    for image, json in dataset2:
        break
    print(json)

def anndata_file_iterator(fileobj, handler=reraise_exception):

    ## KC: would be better if this was streamable.  If it's not possible to read
    ## hdf5/anndata streaming we should either
    ## (1) use a different format (CBOR, TAR, etc) that is streamable
    ## (2) use google primitives to localize (instead of this custom read a stream to disk)
    ##
    filename = "local.h5ad"
    print(f"Localizing File...")
    with open(filename, "bw") as dest:
        while b := fileobj.read(DEFAULT_BUFFER_SIZE):
            dest.write(b)

    adata = anndata.read_h5ad(filename)
    size = adata.raw.X.shape[0]

    for i in range(size):
        try:
            index = adata.obs.index[i]
            data = adata.raw.X[i]
            result = dict(index=index, data=data)
            #print(f"Yielding... {result}")
            yield result
        except Exception as exn:
            if handler(exn):
                continue
            else:
                break
    
    os.remove(filename)


# KC: I find the use of the variable name "sample" here confusing... but this is taken
# from tarexpander as a template
def anndata_file_expander(data, handler=reraise_exception):
    for source in data:
        url = source["url"]
        try:
            assert isinstance(source, dict)
            assert "stream" in source
            for sample in anndata_file_iterator(source["stream"], handler=handler):
                assert (
                    isinstance(sample, dict) and "data" in sample and "index" in sample
                )
                sample["__url__"] = url
                yield sample
        except Exception as exn:
            exn.args = exn.args + (source.get("stream"), source.get("url"))
            if handler(exn):
                continue
            else:
                break

def anndata_samples(src, handler=reraise_exception):
    streams = url_opener(src, handler=handler)
    cells = anndata_file_expander(streams, handler=handler) 
    # we don't need to group by, ids are unique
    # samples = group_by_keys(cells, handler=handler
    return cells


anndata_to_samples = filters.pipelinefilter(anndata_samples)

def use_anndata():
    ## Experiment 2: unroll WebDataset into individual pipeline operations (so we can move to h5ad)
    url = "gs://dsp-cell-annotation-service/benchmark_v1/benchmark_v1.{000..002}.h5ad"

    dataset2 = wds.DataPipeline(
        wds.SimpleShardList(url),
        # at this point we have an iterator over all the shards
        wds.shuffle(100),
        wds.split_by_worker,
        # at this point, we have an iterator over the shards assigned to each worker
        anndata_to_samples(),
        wds.shuffle(1000),
        wds.to_tuple("index", "data")
    )

    isinstance(dataset2, torch.utils.data.IterableDataset)

    i=0
    for x in dataset2:
        #print(x)
        i = i + 1
        #if i == 5:
        #    break
    print(f"processed {i}")

def main():
    #baseline()
    #use_gs_urls()
    #use_data_pipeline()
    use_anndata()

if __name__ == "__main__":
    main()

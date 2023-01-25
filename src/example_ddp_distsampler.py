import torch
from torch import optim, nn, utils, Tensor
import torch.nn as nn
from torch.optim import Adam

from torch.utils.data import DataLoader
from casp.ml.data import DistributedAnnCollection, DistributedAnnCollectionSampler, DistributedAnnCollectionDataset
from lightning.pytorch import Trainer, LightningDataModule

import typing as t


# define the LightningModule
class LitAutoEncoder(LightningModule):
    def __init__(self):
        super().__init__()
        features = 36350
        self.encoder = nn.Sequential(nn.Linear(features, 64), nn.ReLU(), nn.Linear(64, 3))
        self.decoder = nn.Sequential(nn.Linear(3, 64), nn.ReLU(), nn.Linear(64, features))

    def training_step(self, batch, batch_idx):
        # training_step defines the train loop.
        # it is independent of forward
        x, y = batch
        x = x.view(x.size(0), -1)
        z = self.encoder(x)
        x_hat = self.decoder(z)
        loss = nn.functional.mse_loss(x_hat, x)
        # Logging to TensorBoard by default
        self.log("train_loss", loss)
        return loss

    def configure_optimizers(self):
        optimizer = optim.Adam(self.parameters(), lr=1e-3)
        return optimizer

class CustomDataModule(LightningDataModule):
    def __init__(self, url_pattern: str, batch_size: int = 4):
        super().__init__()
        self.url_pattern = url_pattern
        self.batch_size = batch_size

    def train_dataloader(self):

        
        dataset = DistributedAnnCollectionDataset(DistributedAnnCollection(self.url_pattern))

        # debugging without a sampler
        # train_loader = DataLoader(self.dataset, batch_size=self.batch_size)

        # using normal distributedsampler
        # train_sampler = torch.utils.data.distributed.DistributedSampler(dataset, shuffle=True)

        train_sampler = DistributedAnnCollectionSampler(dataset=dataset)

        train_loader = utils.data.DataLoader(
                                        dataset, 
                                        batch_size=self.batch_size,
                                        num_workers=2,
                                        persistent_workers=True,
                                        pin_memory=True,                                        
                                        sampler=train_sampler
                                        )
 
        return train_loader

def main():
    model = LitAutoEncoder()

    batch_size = 100
    url_pattern = "gs://dsp-cell-annotation-service/benchmark_v1/benchmark_v1.{000..100}.h5ad"

    dm = CustomDataModule(url_pattern, batch_size)

    trainer = Trainer(accelerator="gpu", 
                      devices=2, 
                      replace_sampler_ddp=False, 
                      reload_dataloaders_every_n_epochs=1, 
                      max_epochs=10)
    trainer.fit(model, dm)

if __name__ == '__main__':
    main()

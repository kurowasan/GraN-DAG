import os

import numpy as np
import torch


class DataManagerFile(object):
    def __init__(self, file_path, i_dataset, num_folds=1, fold=0, train=True, normalize=False, mean=None, std=None):
        adjacency = np.load(os.path.join(file_path, "DAG{}.npy".format(i_dataset)))
        self.adjacency = torch.as_tensor(adjacency).type(torch.Tensor)
        self.cpdag = np.load(os.path.join(file_path, "CPDAG{}.npy".format(i_dataset)))
        self.mean, self.std = mean, std
        data = np.load(os.path.join(file_path, "data{}.npy".format(i_dataset)))
        total_num_samples = data.shape[0]
        assert 0 <= fold <= num_folds - 1
        fold_size = total_num_samples // num_folds
        if not train:
            self.dataset = torch.as_tensor(data[fold * fold_size : (fold + 1) * fold_size]).type(torch.Tensor)
        else:
            data1, data2 = data[:fold * fold_size], data[(fold + 1) * fold_size: ]
            self.dataset = torch.as_tensor(np.concatenate([data1, data2], 0)).type(torch.Tensor)

        if normalize:
            if self.mean is None or self.std is None:
                self.mean = torch.mean(self.dataset, 0, keepdim=True)
                self.std = torch.std(self.dataset, 0, keepdim=True)
            self.dataset = (self.dataset - self.mean) / self.std

        self.num_samples = self.dataset.size(0)

    def sample(self, batch_size):
        sample_idxs = np.random.choice(np.arange(int(self.num_samples)), size=(int(batch_size),), replace=False)
        samples = self.dataset[torch.as_tensor(sample_idxs).long()]

        return samples, torch.ones_like(samples)  # second output is mask (for intervention in the future)


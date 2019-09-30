# coding=utf-8
"""
GraN-DAG

Copyright © 2019 Authors of Gradient-Based Neural DAG Learning

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated
documentation files (the "Software"), to deal in the Software without restriction, including without limitation the
rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit
persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the
Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE
WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
"""
import os

import numpy as np
import torch


class DataManagerFile(object):
    def __init__(self, file_path, i_dataset, train_samples=0.8, test_samples=None, train=True, normalize=False,
                 mean=None, std=None, random_seed=42):
        """
        Parameters:
        -----------
        train_samples: uint or float, default=0.8
            If float, specifies the proportion of data used for training and the rest is used for testing. If an
            integer, specifies the exact number of examples to use for training.
        test_samples: uint, default=None
            Specifies the number of examples to use for testing. The default value uses all examples that are not used
            for training.

        """
        self.random = np.random.RandomState(random_seed)

        # Load the graph
        adjacency = np.load(os.path.join(file_path, "DAG{}.npy".format(i_dataset)))
        self.adjacency = torch.as_tensor(adjacency).type(torch.Tensor)

        # Load data
        self.data_path = os.path.join(file_path, "data{}.npy".format(i_dataset))
        data = np.load(self.data_path)

        # Determine train/test partitioning
        if isinstance(train_samples, float):
            train_samples = int(data.shape[0] * train_samples)
        if test_samples is None:
            test_samples = data.shape[0] - train_samples
        assert train_samples + test_samples <= data.shape[0], "The number of examples to load must be smaller than " + \
            "the total size of the dataset"

        # Shuffle and filter examples
        shuffle_idx = np.arange(data.shape[0])
        self.random.shuffle(shuffle_idx)
        data = data[shuffle_idx[: train_samples + test_samples]]

        # Train/test split
        if not train:
            if train_samples == data.shape[0]: # i.e. no test set
                self.dataset = None
            else:
                self.dataset = torch.as_tensor(data[train_samples: train_samples + test_samples]).type(torch.Tensor)
        else:
            self.dataset = torch.as_tensor(data[: train_samples]).type(torch.Tensor)

        # Normalize data
        self.mean, self.std = mean, std
        if normalize:
            if self.mean is None or self.std is None:
                self.mean = torch.mean(self.dataset, 0, keepdim=True)
                self.std = torch.std(self.dataset, 0, keepdim=True)
            self.dataset = (self.dataset - self.mean) / self.std

        self.num_samples = self.dataset.size(0)

    def sample(self, batch_size):
        sample_idxs = self.random.choice(np.arange(int(self.num_samples)), size=(int(batch_size),), replace=False)
        samples = self.dataset[torch.as_tensor(sample_idxs).long()]
        return samples, torch.ones_like(samples)  # second output is mask (for intervention in the future)


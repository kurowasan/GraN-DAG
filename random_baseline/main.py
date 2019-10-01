import os
import time
import cdt
import argparse
import pandas as pd
import numpy as np
import networkx as nx
import torch

import sys
sys.path.append("..")

from gran_dag.plot import plot_adjacency
from gran_dag.utils.save import dump
from gran_dag.utils.metrics import edge_errors
from gran_dag.dag_optim import is_acyclic
from gran_dag.data import DataManagerFile
from gran_dag.train import cam_pruning_, pns_


def sample_random_dag(num_nodes, connectProb):
    adj = np.zeros((num_nodes, num_nodes))
    for i in range(1, num_nodes):
        for j in range(i):
            u = np.random.uniform(0, 1)
            adj[i, j] = float((u <= connectProb))

    order = np.random.permutation(num_nodes)
    adj = adj[order, :]
    adj = adj[:, order]
    return adj


def main(opt, metrics_callback=None, plotting_callback=None):
    # Control as much randomness as possible
    torch.manual_seed(opt.random_seed)
    np.random.seed(opt.random_seed)
    torch.set_default_tensor_type('torch.DoubleTensor')

    time0 = time.time()
    opt.model = "random"

    # load data
    train_data = DataManagerFile(opt.data_path, opt.i_dataset, opt.train_samples, opt.test_samples, train=True,
                                 normalize=opt.normalize_data, random_seed=opt.random_seed)
    gt_dag = train_data.adjacency.detach().cpu().numpy()
    train_data_np = train_data.dataset.detach().cpu().numpy()
    num_nodes = train_data_np.shape[1]

    connectProb = np.random.uniform(0, 1)
    adj = sample_random_dag(num_nodes, connectProb)

    # Compute graph metrics
    sid = float(cdt.metrics.SID(target=gt_dag, pred=adj))
    shd = float(cdt.metrics.SHD(target=gt_dag, pred=adj, double_for_anticausal=False))
    shd_cpdag = float(cdt.metrics.SHD_CPDAG(target=gt_dag, pred=adj))
    fn, fp, rev = edge_errors(adj, gt_dag)
    timing = time.time() - time0

    #save
    if not os.path.exists(opt.exp_path):
        os.makedirs(opt.exp_path)

    if metrics_callback is not None:
        metrics_callback(stage="random", step=0,
                         metrics={"sid": sid, "shd": shd,
                                  "shd_cpdag": shd_cpdag, "fn": fn, "fp": fp, "rev": rev},
                         throttle=False)

    dump(opt, opt.exp_path, 'opt')
    dump(timing, opt.exp_path, 'timing', True)
    dump(sid, opt.exp_path, 'sid', True)
    dump(shd, opt.exp_path, 'shd', True)
    np.save(os.path.join(opt.exp_path, "DAG"), adj)

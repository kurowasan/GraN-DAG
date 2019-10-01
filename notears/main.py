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
from notears.notears import notears, retrain

def main(opt, metrics_callback, plotting_callback=None):
    # Control as much randomness as possible
    torch.manual_seed(opt.random_seed)
    np.random.seed(opt.random_seed)
    torch.set_default_tensor_type('torch.DoubleTensor')

    time0 = time.time()
    opt.model = "notears"

    # load data
    train_data = DataManagerFile(opt.data_path, opt.i_dataset, opt.train_samples, opt.test_samples, train=True,
                                 normalize=opt.normalize_data, random_seed=opt.random_seed)
    test_data = DataManagerFile(opt.data_path, opt.i_dataset, opt.train_samples, opt.test_samples, train=False,
                                normalize=opt.normalize_data, mean=train_data.mean, std=train_data.std,
                                random_seed=opt.random_seed)
    gt_dag = train_data.adjacency.detach().cpu().numpy()
    train_data_np = train_data.dataset.detach().cpu().numpy()
    num_nodes = train_data_np.shape[1]

    initial_adj = np.ones((num_nodes, num_nodes)) - np.eye(num_nodes)

    if opt.pns:
        initial_adj = pns_(initial_adj, train_data, test_data, opt.num_neighbors, opt.pns_thresh)

    # apply NOTEARS algorithm
    w_adj, flag_max_iter = notears(train_data_np, initial_adj, opt.lambda1, opt.max_iter, opt.h_tol, opt.w_threshold)

    # make sure is acyclic
    adj = (w_adj != 0).astype(np.double)
    original_adj_cyclic = not is_acyclic(adj)
    while not is_acyclic(adj):
        print("Removing an edge since original DAG was not acyclic")
        w_adj_abs = np.abs(w_adj)
        min_abs_value = np.min(w_adj_abs[np.nonzero(w_adj_abs)])
        to_keep = (w_adj_abs > min_abs_value).astype(np.double)
        w_adj = w_adj * to_keep
        adj = (w_adj != 0).astype(np.double)

    # cam pruning?
    if opt.cam_pruning:
        new_adj = cam_pruning_(adj, train_data, test_data, opt.cutoff, opt.exp_path)
        assert (adj >= new_adj).all()  # assert that cam_pruning is not adding edges
        adj = new_adj

    # evaluate held-out likelihood
    score_train_noretrain = -0.5 * np.sum((train_data_np - np.matmul(train_data_np, adj * w_adj)) ** 2) / train_data_np.shape[0]
    if test_data.dataset is not None:
        test_data_np = test_data.dataset.detach().cpu().numpy()
        score_valid_noretrain = -0.5 * np.sum((test_data_np - np.matmul(test_data_np, adj * w_adj)) ** 2) / test_data_np.shape[0]
    else:
        score_valid_noretrain= False

    # retrain to evaluate held-out likelihood
    score_train, score_valid, flag_max_iter_retrain = retrain(adj, train_data, test_data, opt.max_iter_retrain)

    print("without retrain: train: {}, valid: {}".format(score_train_noretrain, score_valid_noretrain))

    # Compute graph metrics
    sid = float(cdt.metrics.SID(target=gt_dag, pred=adj))
    shd = float(cdt.metrics.SHD(target=gt_dag, pred=adj, double_for_anticausal=False))
    shd_cpdag = float(cdt.metrics.SHD_CPDAG(target=gt_dag, pred=adj))
    fn, fp, rev = edge_errors(adj, gt_dag)
    timing = time.time() - time0

    #save
    if not os.path.exists(opt.exp_path):
        os.makedirs(opt.exp_path)

    metrics_callback(stage="notears", step=0,
                     metrics={"train_score": score_train, "val_score": score_valid, "sid": sid, "shd": shd,
                              "shd_cpdag": shd_cpdag, "fn": fn, "fp": fp, "rev": rev,
                              "original_adj_cyclic": original_adj_cyclic, "flag_max_iter": flag_max_iter,
                              "flag_max_iter_retrain": flag_max_iter_retrain},
                     throttle=False)

    dump(opt, opt.exp_path, 'opt')
    dump(timing, opt.exp_path, 'timing', True)
    dump(score_train, opt.exp_path, 'train_score', True)
    dump(score_valid, opt.exp_path, 'test_score', True)
    dump(sid, opt.exp_path, 'sid', True)
    dump(shd, opt.exp_path, 'shd', True)
    np.save(os.path.join(opt.exp_path, "DAG"), adj)

def _print_metrics(stage, step, metrics, throttle=None):
    for k, v in metrics.items():
        print("    %s:" % k, v)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument('--data-path', type=str, default=None,
                        help='Path to data files')
    parser.add_argument('--i-dataset', type=str, default=None,
                        help='dataset index')
    parser.add_argument('--exp-path', type=str, default='exp',
                        help='Path to experiments')
    parser.add_argument('--pns', action='store_true')
    parser.add_argument('--cam-pruning', action='store_true')
    parser.add_argument('--max-iter-retrain', type=int, default=10000)

    parser.add_argument('--train-samples', type=int, default=0.8,
                        help='Number of samples used for training (default is 80% of the total size)')
    parser.add_argument('--test-samples', type=int, default=None,
                        help='Number of samples used for testing (default is whatever is not used for training)')
    parser.add_argument('--normalize-data', action="store_true",
                        help='(x - mu) / std')
    parser.add_argument('--random-seed', type=int, default=42,
                        help="Random seed for pytorch and numpy")
    parser.add_argument('--lambda1', type=float, default=0.01,
                        help="L1 coeff")
    parser.add_argument('--max-iter', type=int, default=100,
                        help="maximal number of subproblems to solve")
    parser.add_argument('--h-tol', type=float, default=1e-8,
                        help="threshold to declare end of AL")
    parser.add_argument('--w-threshold', type=float, default=0.3,
                        help="Final thresholding")

    opt = parser.parse_args()

    main(opt, metrics_callback=_print_metrics)

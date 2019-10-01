import os
import time
import cdt
import argparse
import pandas as pd
import numpy as np
import networkx as nx

import sys
sys.path.append("..")
from gran_dag.plot import plot_adjacency
from gran_dag.utils.save import dump
from gran_dag.utils.metrics import edge_errors
from gran_dag.data import DataManagerFile
from cam import CAM_with_score


def main(opt, metrics_callback=None, plotting_callback=None):
    time0 = time.time()
    opt.model = "cam"

    # load data
    train_data = DataManagerFile(opt.data_path, opt.i_dataset, opt.train_samples, opt.test_samples, train=True,
                                 normalize=opt.normalize_data, random_seed=opt.random_seed)
    test_data = DataManagerFile(opt.data_path, opt.i_dataset, opt.train_samples, opt.test_samples, train=False,
                                normalize=opt.normalize_data, mean=train_data.mean, std=train_data.std,
                                random_seed=opt.random_seed)
    gt_dag = train_data.adjacency.detach().cpu().numpy()
    train_data_pd = pd.DataFrame(train_data.dataset.detach().cpu().numpy())

    # apply CAM
    if test_data.dataset is not None:
        test_data_pd = pd.DataFrame(test_data.dataset.detach().cpu().numpy())
        obj = CAM_with_score(opt.score, opt.cutoff, opt.variable_sel, opt.sel_method,
                             opt.pruning, opt.prune_method)
        dag, train_score, val_score = obj.get_score(train_data_pd, test_data_pd)
    else:
        obj = cdt.causality.graph.CAM(opt.score, opt.cutoff, opt.variable_sel, opt.sel_method, opt.pruning,
                                      opt.prune_method)
        dag = obj.predict(train_data_pd)
        train_score, val_score = None, None
    dag = nx.to_numpy_matrix(dag)

    # Compute SHD and SID metrics
    sid = float(cdt.metrics.SID(target=gt_dag, pred=dag))
    shd = float(cdt.metrics.SHD(target=gt_dag, pred=dag, double_for_anticausal=False))
    shd_cpdag = float(cdt.metrics.SHD_CPDAG(target=gt_dag, pred=dag))
    fn, fp, rev = edge_errors(dag, gt_dag)
    timing = time.time() - time0

    #save
    if not os.path.exists(opt.exp_path):
        os.makedirs(opt.exp_path)

    if metrics_callback is not None:
        metrics_callback(stage="cam", step=0,
                         metrics={"train_score": train_score, "val_score": val_score, "sid": sid, "shd": shd,
                                  "shd_cpdag": shd_cpdag, "fn": fn, "fp": fp, "rev": rev}, throttle=False)

    dump(opt, opt.exp_path, 'opt')
    dump(timing, opt.exp_path, 'timing', True)
    dump(train_score, opt.exp_path, 'train_score', True)
    dump(val_score, opt.exp_path, 'test_score', True)
    dump(sid, opt.exp_path, 'sid', True)
    dump(shd, opt.exp_path, 'shd', True)
    np.save(os.path.join(opt.exp_path, "DAG"), dag)

    plot_adjacency(gt_dag, dag, opt.exp_path)

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
    # parser.add_argument('--score', type=str, default='nonlinear',
    #                     help='Score used to fit the gaussian processes')

    # Variable selection (PNS)
    parser.add_argument('--variable-sel', action="store_true",
                        help='Perform a variable selection step')
    # parser.add_argument('--sel-method', type=str, default='gamboost',
    #                     help='Method used for variable selection')
    parser.add_argument('--cutoff', type=float, default=0.001,
                        help='Threshold value for vaiable selection')

    # Pruning
    parser.add_argument('--pruning', action="store_true", default=True,
                        help='Perform an initial pruning step')
    # parser.add_argument('--prune-method', type=str, default="gam",
    #                     help='after to-dag or pruning, retrain model from scratch before reporting nll-val')
    parser.add_argument('--train-samples', type=int, default=0.8,
                        help='Number of samples used for training (default is 80% of the total size)')
    parser.add_argument('--test-samples', type=int, default=None,
                        help='Number of samples used for testing (default is whatever is not used for training)')
    parser.add_argument('--normalize-data', action="store_true",
                        help='(x - mu) / std')
    parser.add_argument('--random-seed', type=int, default=42,
                        help="Random seed for pytorch and numpy")

    opt = parser.parse_args()
    opt.score = 'nonlinear'
    opt.sel_method = 'gamboost'
    opt.prune_method = 'gam'

    main(opt, metrics_callback=_print_metrics)

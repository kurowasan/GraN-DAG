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

def main(opt, metrics_callback, plotting_callback=None):
    time0 = time.time()
    opt.model = "pc"

    # load data
    train_data = DataManagerFile(opt.data_path, opt.i_dataset, opt.train_samples, opt.test_samples, train=True,
                                 normalize=opt.normalize_data, random_seed=opt.random_seed)
    test_data = DataManagerFile(opt.data_path, opt.i_dataset, opt.train_samples, opt.test_samples, train=False,
                                normalize=opt.normalize_data, mean=train_data.mean, std=train_data.std,
                                random_seed=opt.random_seed)
    gt_dag = train_data.adjacency.detach().cpu().numpy()
    if test_data is not None:
        all_data_np = np.concatenate([train_data.dataset.detach().cpu().numpy(), test_data.dataset.detach().cpu().numpy()], 0)
    else:
        all_data_np = train_data.dataset.detach().cpu().numpy()
    all_data_pd = pd.DataFrame(all_data_np)

    # apply PC
    obj = cdt.causality.graph.PC(opt.ci_test, opt.method_indep, opt.alpha)
    dag = obj.predict(all_data_pd)
    train_score, val_score = None, None
    dag = nx.to_numpy_matrix(dag)

    # Compute SHD and SID metrics
    sid = float(cdt.metrics.SID(target=gt_dag, pred=dag))
    shd = None
    shd_cpdag = float(cdt.metrics.SHD_CPDAG(target=gt_dag, pred=dag))
    fn, fp, rev = None, None, None
    timing = time.time() - time0

    #save
    if not os.path.exists(opt.exp_path):
        os.makedirs(opt.exp_path)

    metrics_callback(stage="pc", step=0,
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


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument('--data-path', type=str, default=None,
                        help='Path to data files')
    parser.add_argument('--i-dataset', type=str, default=None,
                        help='dataset index')
    parser.add_argument('--exp-path', type=str, default='exp',
                        help='Path to experiments')

    # PC parameters
    parser.add_argument('--ci-test', type=str, default='gaussian',
                        help='Conditional independence test (gaussian| \
                        hsic|discrete|binary|randomized)')
    parser.add_argument('--method-indep', type=str, default='corr',
                        help='Heuristic for testing CI test (dcc| \
                        hsic_gamma|hsic_perm|hsic_clust|corr|rcit|rcot)')
    parser.add_argument('--alpha', type=float, default=0.01,
                        help='significance level for the individual ci test')

    opt = parser.parse_args()
    main(opt)

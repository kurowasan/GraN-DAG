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
from ges import GES_with_score


def main(opt, metrics_callback, plotting_callback=None):
    time0 = time.time()
    opt.model = "ges"

    # load data
    # TODO: Use same dataloader as gran_dag to make sure we have the same train/valid split
    # load data
    train_data = DataManagerFile(opt.data_path, opt.i_dataset, opt.train_samples, opt.test_samples, train=True,
                                 normalize=opt.normalize_data, random_seed=opt.random_seed)
    test_data = DataManagerFile(opt.data_path, opt.i_dataset, opt.train_samples, opt.test_samples, train=False,
                                normalize=opt.normalize_data, mean=train_data.mean, std=train_data.std,
                                random_seed=opt.random_seed)
    gt_dag = train_data.adjacency.detach().cpu().numpy()
    train_data_pd = pd.DataFrame(train_data.dataset.detach().cpu().numpy())

    # apply GES
    if test_data.dataset is not None:
        test_data_pd = pd.DataFrame(test_data.dataset.detach().cpu().numpy())
        obj = GES_with_score(opt.lambda_ges)
        dag, train_score, val_score = obj.get_score(train_data_pd, test_data_pd)
    else:
        obj = cdt.causality.graph.GES()
        dag = obj.predict(train_data_pd)
        train_score, val_score = None, None
    dag = nx.to_numpy_matrix(dag)

    # Compute SHD-CPDAG and SID metrics
    sid = float(cdt.metrics.SID(target=gt_dag, pred=dag))
    shd = None # Since GES outputs a CPDAG
    shd_cpdag = float(cdt.metrics.SHD_CPDAG(target=gt_dag, pred=dag))
    fn, fp, rev = None, None, None
    timing = time.time() - time0

    #save
    if not os.path.exists(opt.exp_path):
        os.makedirs(opt.exp_path)

    metrics_callback(stage="ges", step=0,
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
    # lambda should be equal to log(n)/2 where n is the number of examples
    parser.add_argument('--lambda-ges', type=float, default=1,
                        help='Penalization constant used by GES')

    opt = parser.parse_args()
    main(opt)

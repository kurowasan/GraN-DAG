import networkx as nx
import numpy as np
import os
import warnings
import multiprocessing as mp
from functools import partial

import cppext

from data import load_from_folder
from notears import notears
from tqdm import tqdm

def hyperparameter_search(lambdas_1, w_thresholds, data, num_cpus):
    if not np.all(np.diff(w_thresholds) > 0):
        raise ValueError('The w_thresholds must be in increasing order')
    # Use the first experiment to do the hyperparameter search
    graph_ground_truth = nx.DiGraph(data.dag)
    d = len(graph_ground_truth)

    with mp.Pool(num_cpus) as pool:
        results = pool.map(partial(notears, data.train, max_iter=args.max_iter,
            h_tol=args.h_tol, w_threshold=args.w_threshold, G=graph_ground_truth), lambdas_1)

    valid_losses = np.zeros((len(lambdas_1), len(w_thresholds)))
    is_dag = np.zeros((len(lambdas_1), len(w_thresholds)), dtype=np.bool_)
    for i, lambda1 in enumerate(lambdas_1):
        w_est_sparse = np.copy(results[i].w_est).flatten()
        for j, w_threshold in enumerate(w_thresholds):
            w_est_sparse[np.abs(w_est_sparse) < w_threshold] = 0
            valid_losses[i, j] = cppext.F_func(w_est_sparse, data.valid, lambda1)

            graph_estimated = nx.DiGraph(w_est_sparse.reshape((d, d)))
            is_dag[i, j] = nx.is_directed_acyclic_graph(graph_estimated)

    x, y = np.nonzero(is_dag)
    index = np.argmin(valid_losses[is_dag])
    lambda1, w_threshold = lambdas_1[x[index]], w_thresholds[y[index]]

    with open(os.path.join(args.output_folder, 'hyperparams.npz'), 'wb') as f:
        np.savez(f, valid_losses=valid_losses, is_dag=is_dag, lambdas_1=lambdas_1,
            w_thresholds=w_thresholds, lambda1=lambda1, w_threshold=w_threshold)

    return lambda1, w_threshold

def no_tears_for_data(args, _args):
    i, data = _args

    graph_ground_truth = nx.DiGraph(data.dag)
    if not nx.is_directed_acyclic_graph(graph_ground_truth):
        warnings.warn('The ground truth graph ({0}) is not acyclic. Ignoring.'.format(i))
        return False

    result = notears(data.train, args.lambda1, max_iter=args.max_iter,
        h_tol=args.h_tol, w_threshold=args.w_threshold, G=graph_ground_truth)

    with open(os.path.join(args.output_folder, 'DAG{0}.npy'.format(i)), 'wb') as f:
        np.save(f, result.w_est_sparse)

    with open(os.path.join(args.output_folder, 'RES{0}.npz'.format(i)), 'wb') as f:
        graph_estimated = nx.DiGraph(result.w_est_sparse)
        np.savez(f, lambda1=args.lambda1, max_iter=args.max_iter, h_tol=args.h_tol,
            w_threshold=args.w_threshold, F_true=result.F_true, w_est=result.w_est,
            w_est_sparse=result.w_est_sparse, iters=result.iters, F_est=result.F_est,
            is_dag=int(nx.is_directed_acyclic_graph(graph_estimated)))

    return True

def main(args):
    datas = load_from_folder(args.folder, valid_split=args.valid_split)
    num_cpus = int(os.environ.get('SLURM_CPUS_PER_TASK', 1))
    if args.num_experiments < 0:
        args.num_experiments = len(datas)

    if args.hyperparams_search:
        lambda1, w_threshold = hyperparameter_search([0, 1e-4, 5e-4, 1e-3, 5e-3, 1e-2, 5e-2, 0.1],
            [0.003, 0.01, 0.03, 0.1, 0.3], datas[0], num_cpus)
        args.lambda1 = lambda1
        args.w_threshold = w_threshold

    with mp.Pool(num_cpus) as pool:
        list(tqdm(pool.imap(partial(no_tears_for_data, args), enumerate(datas[:args.num_experiments])),
            total=args.num_experiments))

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser('DAGS with NOTEARS')

    parser.add_argument('folder', type=str,
        help='Path to the data folder')
    parser.add_argument('--output-folder', type=str,
        help='Path to the output folder')
    parser.add_argument('--lambda1', type=float, default=0.01,
        help='Value of the l1 regularization parameter')
    parser.add_argument('--valid-split', type=float, default=0.2,
        help='Amount of data for validation')
    parser.add_argument('--max_iter', type=int, default=100,
        help='Maximum number of dual ascent steps')
    parser.add_argument('--h_tol', type=float, default=1e-8,
        help='Exit if |h(w)| <= h_tol')
    parser.add_argument('--w_threshold', type=float, default=0.3,
        help='Fixed threshold for edge weights')
    parser.add_argument('--hyperparams-search', action='store_true',
        help='Do an hyperparameter search')
    parser.add_argument('--num-experiments', type=int, default=-1,
        help='Number of experiments')

    args = parser.parse_args()
    main(args)

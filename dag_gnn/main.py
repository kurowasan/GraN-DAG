import sys
import os
import time
import cdt
import argparse
import pandas as pd
import numpy as np
import torch
import networkx as nx

sys.path.append("..")
from gran_dag.plot import plot_adjacency
from gran_dag.utils.save import dump
from gran_dag.utils.metrics import edge_errors
from gran_dag.dag_optim import is_acyclic
from gran_dag.data import DataManagerFile
from gran_dag.train import cam_pruning_, pns_
from dag_gnn.train import dag_gnn, retrain
from dag_gnn.utils import load_numpy_data

def main(opt, metrics_callback, plotting_callback=None):
    # Control as much randomness as possible
    torch.manual_seed(opt.random_seed)
    np.random.seed(opt.random_seed)
    torch.set_default_tensor_type('torch.DoubleTensor')

    time0 = time.time()
    opt.model = "dag_gnn"

    # load data
    train_data = DataManagerFile(opt.data_path, opt.i_dataset, opt.train_samples, opt.test_samples, train=True,
                                 normalize=opt.normalize_data, random_seed=opt.random_seed)
    test_data = DataManagerFile(opt.data_path, opt.i_dataset, opt.train_samples, opt.test_samples, train=False,
                                normalize=opt.normalize_data, mean=train_data.mean, std=train_data.std,
                                random_seed=opt.random_seed)
    gt_dag = train_data.adjacency.detach().cpu().numpy()
    train_data_np = train_data.dataset.unsqueeze(2).detach().cpu().numpy()  # shape (num_samples, num_nodes, node_dim)
    num_nodes = train_data_np.shape[1]

    initial_adj = np.ones((num_nodes, num_nodes), dtype=np.double) - np.eye(num_nodes, dtype=np.double)

    if opt.pns:
        initial_adj = pns_(initial_adj, train_data, test_data, opt.num_neighbors, opt.pns_thresh)

    # apply DAG-GNN
    w_adj, rel_rec, rel_send, encoder, decoder, train_loader, flag_max_iter = dag_gnn(opt, train_data_np, gt_dag, initial_adj)
    adj = (w_adj != 0).astype(np.double)

    # verify the masking worked properly
    assert (initial_adj >= adj).all() # assert that edge_init = 0 => edge = 0

    # make sure is acyclic
    original_adj_cyclic = not is_acyclic(adj)
    while not is_acyclic(adj):
        print("Removing an edge since original DAG was not acyclic")
        w_adj_abs = np.abs(w_adj)
        min_abs_value = np.min(w_adj_abs[np.nonzero(w_adj_abs)])
        to_keep = (w_adj_abs > min_abs_value).astype(np.double)
        w_adj = w_adj * to_keep
        adj = (w_adj != 0).astype(np.double)

    # cam_pruning?
    if opt.cam_pruning:
        new_adj = cam_pruning_(adj, train_data, test_data, opt.cutoff, opt.exp_path)
        assert (adj >= new_adj).all() # assert that cam_pruning is not adding edges
        adj = new_adj

    # evaluate held-out likelihood
    if test_data.dataset is not None:
        test_data_np = test_data.dataset.unsqueeze(2).detach().cpu().numpy()
    else:
        test_data_np = None
    score_train, score_valid, flag_max_iter_retrain = retrain(opt, train_data_np, test_data_np, gt_dag, adj)

    # Compute SHD and SID metrics
    sid = float(cdt.metrics.SID(target=gt_dag, pred=adj))
    shd = float(cdt.metrics.SHD(target=gt_dag, pred=adj, double_for_anticausal=False))
    shd_cpdag = float(cdt.metrics.SHD_CPDAG(target=gt_dag, pred=adj))
    fn, fp, rev = edge_errors(adj, gt_dag)
    timing = time.time() - time0

    #save
    if not os.path.exists(opt.exp_path):
        os.makedirs(opt.exp_path)

    metrics_callback(stage="dag_gnn", step=0,
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

    #plot_adjacency(gt_dag, adj, opt.exp_path)

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
    parser.add_argument('--pns-thres', type=float, default=0.75)
    parser.add_argument('--cutoff', type=float, default=0.001)

    # -----------data parameters ------
    parser.add_argument('--data_sample_size', type=int, default=1000,
                        help='the number of samples of data')
    parser.add_argument('--data_variable_size', type=int, default=10,
                        help='the number of variables in synthetic generated data')
    parser.add_argument('--x_dims', type=int, default=1,  # changed here
                        help='The number of input dimensions: default 1.')
    parser.add_argument('--z_dims', type=int, default=1,
                        help='The number of latent variable dimensions: default the same as variable size.')

    # -----------training hyperparameters
    parser.add_argument('--optimizer', type=str, default='Adam',
                        help='the choice of optimizer used')
    parser.add_argument('--graph_threshold', type=float, default=0.3,  # 0.3 is good, 0.2 is error prune
                        help='threshold for learned adjacency matrix binarization')
    parser.add_argument('--tau_A', type=float, default=0.0,
                        help='coefficient for L-1 norm of A.')
    parser.add_argument('--lambda_A', type=float, default=0.,
                        help='coefficient for DAG constraint h(A).')
    parser.add_argument('--c_A', type=float, default=1,
                        help='coefficient for absolute value h(A).')
    parser.add_argument('--use_A_connect_loss', type=int, default=0,
                        help='flag to use A connect loss')
    parser.add_argument('--use_A_positiver_loss', type=int, default=0,
                        help='flag to enforce A must have positive values')

    parser.add_argument('--no-cuda', action='store_true', default=True,
                        help='Disables CUDA training.')
    #parser.add_argument('--seed', type=int, default=42, help='Random seed.') not used anywhere...
    parser.add_argument('--epochs', type=int, default=200,
                        help='Number of epochs to train.')
    parser.add_argument('--batch-size', type=int, default=100,
                        # note: should be divisible by sample size, otherwise throw an error
                        help='Number of samples per batch.')
    parser.add_argument('--lr', type=float, default=3e-3,  # basline rate = 1e-3
                        help='Initial learning rate.')
    parser.add_argument('--encoder-hidden', type=int, default=64,
                        help='Number of hidden units.')
    parser.add_argument('--decoder-hidden', type=int, default=64,
                        help='Number of hidden units.')
    parser.add_argument('--temp', type=float, default=0.5,
                        help='Temperature for Gumbel softmax.')
    parser.add_argument('--k_max_iter', type=int, default=1e2,
                        help='the max iteration number for searching lambda and c')

    parser.add_argument('--encoder', type=str, default='mlp',
                        help='Type of path encoder model (mlp, or sem).')
    parser.add_argument('--decoder', type=str, default='mlp',
                        help='Type of decoder model (mlp, or sim).')
    parser.add_argument('--no-factor', action='store_true', default=False,
                        help='Disables factor graph model.')
    parser.add_argument('--suffix', type=str, default='_springs5',
                        help='Suffix for training data (e.g. "_charged".')
    parser.add_argument('--encoder-dropout', type=float, default=0.0,
                        help='Dropout rate (1 - keep probability).')
    parser.add_argument('--decoder-dropout', type=float, default=0.0,
                        help='Dropout rate (1 - keep probability).')
    parser.add_argument('--save-folder', type=str, default='logs',
                        help='Where to save the trained model, leave empty to not save anything.')
    parser.add_argument('--load-folder', type=str, default='',
                        help='Where to load the trained model if finetunning. ' +
                             'Leave empty to train from scratch')

    parser.add_argument('--h_tol', type=float, default=1e-8,
                        help='the tolerance of error of h(A) to zero')
    parser.add_argument('--prediction-steps', type=int, default=10, metavar='N',
                        help='Num steps to predict before re-using teacher forcing.')
    parser.add_argument('--lr-decay', type=int, default=200,
                        help='After how epochs to decay LR by a factor of gamma.')
    parser.add_argument('--gamma', type=float, default=1.0,
                        help='LR decay factor.')
    parser.add_argument('--skip-first', action='store_true', default=False,
                        help='Skip first edge type in decoder, i.e. it represents no-edge.')
    parser.add_argument('--var', type=float, default=5e-5,
                        help='Output variance.')
    parser.add_argument('--hard', action='store_true', default=False,
                        help='Uses discrete samples in training forward pass.')
    parser.add_argument('--prior', action='store_true', default=False,
                        help='Whether to use sparsity prior.')
    parser.add_argument('--dynamic-graph', action='store_true', default=False,
                        help='Whether test with dynamically re-computed graph.')

    parser.add_argument('--train-samples', type=int, default=0.8,
                        help='Number of samples used for training (default is 80% of the total size)')
    parser.add_argument('--test-samples', type=int, default=None,
                        help='Number of samples used for testing (default is whatever is not used for training)')
    parser.add_argument('--normalize-data', action="store_true",
                        help='(x - mu) / std')
    parser.add_argument('--random-seed', type=int, default=42,
                        help="Random seed for pytorch and numpy")
    parser.add_argument('--cam-pruning', action='store_true')


    opt = parser.parse_args()

    main(opt, metrics_callback=_print_metrics)

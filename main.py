import os
import argparse
import copy

from comet_ml import Experiment, ExistingExperiment
import numpy as np
import torch

from models.learnables import LearnableModel_NonLinGauss, LearnableModel_NonLinGaussANM
from train import pns, train, to_dag
from data import DataManagerFile
from utils import load, dump


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    # experiment
    parser.add_argument('--exp-path', type=str, default='/exp',
                        help='Path to experiments')
    parser.add_argument('--pns', action="store_true",
                        help='Run `pns` function, get /pns folder')
    parser.add_argument('--train', action="store_true",
                        help='Run `train` function, get /train folder')
    parser.add_argument('--to-dag', action="store_true",
                        help='Run `to_dag` function, get /to-dag folder')

    # data
    parser.add_argument('--data-path', type=str, default=None,
                        help='Path to data files')
    parser.add_argument('--i-dataset', type=str, default=None,
                        help='dataset index')
    parser.add_argument('--num-vars', type=int, default=2,
                        help='Number of variables')
    parser.add_argument('--num-folds', type=int, default=5,
                        help='number of folds for cross-validation')
    parser.add_argument('--fold', type=int, default=0,
                        help='fold we should use for testing')
    parser.add_argument('--train-batch-size', type=int, default=64,
                        help='number of samples in a minibatch')
    parser.add_argument('--num-train-iter', type=int, default=100000,
                        help='number of meta gradient steps')
    parser.add_argument('--normalize-data', action="store_true",
                        help='(x - mu) / std')

    # model
    parser.add_argument('--model', type=str, required=True,
                        help='model class')
    parser.add_argument('--num-layers', type=int, default=2,
                        help="number of hidden layers")
    parser.add_argument('--hid-dim', type=int, default=10,
                        help="number of hidden units per layer")
    parser.add_argument('--nonlin', type=str, default='leaky-relu',
                        help="leaky-relu | sigmoid")
    parser.add_argument('--num-neighbors', type=int, default=None,
                        help='number of neighbors to select in PNS')
    parser.add_argument('--pns-thresh', type=float, default=0.75,
                        help='threshold in PNS')
    parser.add_argument('--edge-clamp-range', type=float, default=1e-4,
                        help='clamping the edges (i,j) to zero when prod_ij is that close to zero')

    # optimization
    parser.add_argument('--lr', type=float, default=1e-2,
                        help='learning rate for optim')
    parser.add_argument('--lr-reinit', type=float, default=1e-4,
                        help='learning rate for optim after first subproblem')
    parser.add_argument('--stop-crit-win', type=int, default=100,
                        help='window size to compute stopping criterion')

    # Augmented Lagrangian options
    parser.add_argument('--omega-lambda', type=float, default=1e-4,
                        help='Precision to declare convergence of subproblems')
    parser.add_argument('--omega-mu', type=float, default=0.9,
                        help='After subproblem solved, h should have reduced by this ratio')
    parser.add_argument('--mu-init', type=float, default=1e-3,
                        help='initial value of mu')
    parser.add_argument('--lambda-init', type=float, default=0.,
                        help='initial value of lambda')
    parser.add_argument('--h-threshold', type=float, default=1e-8,
                        help='stop when h is this close to zero')
    parser.add_argument('--norm-prod', type=str, default="paths",
                        help='how to normalize prod: paths|none')
    parser.add_argument('--jac-thresh', action="store_true",
                        help='threshold using the Jacobian instead of prod')

    # logging
    parser.add_argument('--plot-freq', type=int, default=10000,
                        help='plotting frequency')

    # device and numerical precision
    parser.add_argument('--gpu', action="store_true",
                        help="Use GPU")
    parser.add_argument('--float', action="store_true",
                        help="Use Float precision")

    opt = parser.parse_args()

    # Use GPU
    if opt.gpu:
        if opt.float:
            torch.set_default_tensor_type('torch.cuda.FloatTensor')
        else:
            torch.set_default_tensor_type('torch.cuda.DoubleTensor')
    else:
        if opt.float:
            torch.set_default_tensor_type('torch.FloatTensor')
        else:
            torch.set_default_tensor_type('torch.DoubleTensor')

    # create experiment path
    if not os.path.exists(opt.exp_path):
        os.makedirs(opt.exp_path)

    # create learning model and ground truth model
    if opt.model == "NonLinGauss":
        model = LearnableModel_NonLinGauss(opt.num_vars, opt.num_layers, opt.hid_dim, nonlin=opt.nonlin)
    elif opt.model == "NonLinGaussANM":
        model = LearnableModel_NonLinGaussANM(opt.num_vars, opt.num_layers, opt.hid_dim, nonlin=opt.nonlin)
    else:
        raise ValueError("opt.model has to be in {NonLinGauss, NonLinGaussANM}")

    # create DataManager for training
    train_data = DataManagerFile(opt.data_path, opt.i_dataset, opt.num_folds, opt.fold, train=True,
                                 normalize=opt.normalize_data)
    test_data = DataManagerFile(opt.data_path, opt.i_dataset, opt.num_folds, opt.fold, train=False,
                                normalize=opt.normalize_data, mean=train_data.mean, std=train_data.std)

    # save gt adjacency
    dump(train_data.adjacency.detach().cpu().numpy(), opt.exp_path, 'gt-adjacency')

    # apply preliminary neighborhood selection
    if opt.pns:
        if opt.num_neighbors is None: num_neighbors = opt.num_vars
        else: num_neighbors = opt.num_neighbors
        model = pns(model, train_data, test_data, num_neighbors, opt.pns_thresh, opt.exp_path)

    # train until constraint is sufficiently close to being satisfied
    if opt.train:
        if os.path.exists(os.path.join(opt.exp_path, "pns")):
            model = load(os.path.join(opt.exp_path, "pns"), "model.pkl")
        model = train(model, train_data.adjacency.detach().cpu().numpy(), train_data, test_data, opt)

    # remove edges until we have a DAG
    if opt.to_dag:
        # load model
        assert os.path.exists(os.path.join(opt.exp_path, "train")), \
            "The /train folder is required to run --to_dag. Add --train to the command line"
        model = load(os.path.join(opt.exp_path, "train"), "model.pkl")

        # run
        model = to_dag(model, train_data, test_data, opt)


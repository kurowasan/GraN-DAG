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
import argparse
import copy
import numpy as np

from gran_dag.main import main


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
    parser.add_argument('--cam-pruning', action="store_true",
                        help='Run `cam_pruning` function, get /cam-pruning folder')
    parser.add_argument('--retrain', action="store_true",
                        help='after to-dag or pruning, retrain model from scratch before reporting nll-val')
    parser.add_argument('--random-seed', type=int, default=42, help="Random seed for pytorch and numpy")

    # data
    parser.add_argument('--data-path', type=str, default=None,
                        help='Path to data files')
    parser.add_argument('--i-dataset', type=str, default=None,
                        help='dataset index')
    parser.add_argument('--num-vars', type=int, default=2,
                        help='Number of variables')
    parser.add_argument('--train-samples', type=int, default=0.8,
                        help='Number of samples used for training (default is 80% of the total size)')
    parser.add_argument('--test-samples', type=int, default=None,
                        help='Number of samples used for testing (default is whatever is not used for training)')
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

    # optimization
    parser.add_argument('--optimizer', type=str, default="rmsprop",
                        help='sgd|rmsprop')
    parser.add_argument('--lr', type=float, default=1e-3,
                        help='learning rate for optim')
    parser.add_argument('--lr-reinit', type=float, default=None,
                        help='Learning rate for optim after first subproblem. Default mode reuses --lr.')
    parser.add_argument('--scale-lr-with-mu', action="store_true",
                        help='Scale the learning rate wrt mu in the augmented lagrangian.')
    parser.add_argument('--stop-crit-win', type=int, default=100,
                        help='window size to compute stopping criterion')

    # pns, pruning and thresholding
    parser.add_argument('--pns-thresh', type=float, default=0.75,
                        help='threshold in PNS')
    parser.add_argument('--num-neighbors', type=int, default=None,
                        help='number of neighbors to select in PNS')
    parser.add_argument('--edge-clamp-range', type=float, default=1e-4,
                        help='as we train, clamping the edges (i,j) to zero when prod_ij is that close to zero. '
                             '0 means no clamping. Uses masks on inputs. Once an edge is clamped, no way back.')
    parser.add_argument('--cam-pruning-cutoff', nargs='+',
                        default=np.logspace(-6, 0, 10),
                        help='list of cutoff values. Higher means more edges')

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
                        help='Stop when |h|<X. Zero means stop AL procedure only when h==0. Should use --to-dag even '
                             'with --h-threshold 0 since we might have h==0 while having cycles (due to numerical issues).')

    # misc
    parser.add_argument('--norm-prod', type=str, default="paths",
                        help='how to normalize prod: paths|none')
    parser.add_argument('--square-prod', action="store_true",
                        help="square weights instead of absolute value in prod")
    parser.add_argument('--jac-thresh', action="store_true",
                        help='threshold using the Jacobian instead of prod')
    parser.add_argument('--patience', type=int, default=10,
                        help='Early stopping patience in --retrain.')

    # logging
    parser.add_argument('--plot-freq', type=int, default=10000,
                        help='plotting frequency')
    parser.add_argument('--no-w-adjs-log', action="store_true",
                        help='do not log weighted adjacency (to save RAM). One plot will be missing (A_\phi plot)')

    # device and numerical precision
    parser.add_argument('--gpu', action="store_true",
                        help="Use GPU")
    parser.add_argument('--float', action="store_true",
                        help="Use Float precision")

    main(parser.parse_args())


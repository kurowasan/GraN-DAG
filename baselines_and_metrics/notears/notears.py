"""The full NOTEARS algorithm with l1 regularization.

Defines two functions:
- notears(): the vanilla python implementation of the algorithm

The full NOTEARS algorithm works with both large n and small n.

Code taken from: https://github.com/xunzheng/notears
"""
import numpy as np
import networkx as nx

import cppext

from collections import namedtuple

NotearsData = namedtuple('NotearsData', 'F_true F_est w_est w_est_sparse iters '
    'lambda1 max_iter h_tol w_threshold')


def notears(X, lambda1, max_iter=100, h_tol=1e-8, w_threshold=0.3, G=None):
    """Solve min_W F(W; X) s.t. h(W) = 0 using augmented Lagrangian.

    Args:
        X: [n,d] sample matrix
        lambda1: l1 regularization parameter
        max_iter: max number of dual ascent steps
        h_tol: exit if |h(w)| <= h_tol
        w_threshold: fixed threshold for edge weights
        G: nx.DiGraph ground-truth graph

    Returns:
        W_est: [d,d] estimate
    """
    if np.isfortran(X):
        X = np.ascontiguousarray(X)

    if G is not None:
        w_true = nx.to_numpy_array(G).flatten()
        F_true = cppext.F_func(w_true, X, lambda1)
    else:
        F_true = None

    n, d = X.shape
    w_est, w_new = np.zeros(d * d), np.zeros(d * d)
    rho, alpha, h, h_new = 1.0, 0.0, np.inf, np.inf
    for iters in range(max_iter):
        while rho < 1e+20:
            w_new = cppext.minimize_subproblem(w_est, X, rho, alpha, lambda1)
            h_new = cppext.h_func(w_new)
            if h_new > 0.25 * h:
                rho *= 10
            else:
                break
        w_est, h = w_new, h_new
        alpha += rho * h
        if h <= h_tol:
            break

    w_est_sparse = np.copy(w_est)
    w_est_sparse[np.abs(w_est_sparse) < w_threshold] = 0

    F_est = cppext.F_func(w_est_sparse.flatten(), X, lambda1)

    return NotearsData(F_true=F_true, F_est=F_est, w_est=w_est.reshape((d, d)),
        w_est_sparse=w_est_sparse.reshape((d, d)), iters=iters,
        lambda1=lambda1, max_iter=max_iter, h_tol=h_tol, w_threshold=w_threshold)

import os

import matplotlib
# To avoid displaying the figures
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage.filters import gaussian_filter1d
import torch

def plot_weighted_adjacency(weighted_adjacency, gt_adjacency, exp_path, name="abs-weight-product", mus=None,
                            lambdas=None, iter=None):
    """iter is useful to deal with jacobian, it will interpolate."""
    num_vars = weighted_adjacency.shape[1]
    max_value = 0
    fig, ax1 = plt.subplots()

    for i in range(num_vars):
        for j in range(num_vars):
            if gt_adjacency[i, j]:
                continue
            else:
                color = 'r'
            y = weighted_adjacency[:, i, j]
            num_iter = len(y)
            if iter is not None and len(y) > 1:
                num_iter = iter + 1
                y = np.interp(np.arange(iter + 1), np.linspace(0, iter, num=len(y), dtype=int), y)
            ax1.plot(range(num_iter), y, color, linewidth=1)
            if len(y) > 0: max_value = max(max_value, np.max(y))

    for i in range(num_vars):
        for j in range(num_vars):
            if gt_adjacency[i, j]:
                color = 'g'
            else:
                continue
            y = weighted_adjacency[:, i, j]
            num_iter = len(y)
            #plt.plot(range(len(weighted_adjacency[:, 0, 0])), y, color, alpha=0.1, linewidth=1)
            if iter is not None and len(y) > 1:
                num_iter = iter + 1
                y = np.interp(np.arange(iter + 1), np.linspace(0, iter, num=len(y), dtype=int), y)
            ax1.plot(range(num_iter), y, color, linewidth=1)
            if len(y) > 0: max_value = max(max_value, np.max(y))

    ax1.set_xlabel("Iterations")
    ax1.set_ylabel(name)
    ax1.set_yscale("log")

    if mus is not None or lambdas is not None:
        ax2 = ax1.twinx()
        ax2.set_ylabel(r'$\frac{\mu}{2}$ and $\lambda$', color='blue')
        if mus is not None:
            ax2.plot(range(len(mus)), 0.5 * np.array(mus), color='blue', linestyle="dashed", linewidth=1,
                     label=r"$\frac{\mu}{2}$")
        if lambdas is not None:
            ax2.plot(range(len(lambdas)), lambdas, color='blue', linestyle="dotted", linewidth=1, label=r"$\lambda$")
        ax2.legend()
        ax2.set_yscale("log")
        ax2.tick_params(axis='y', labelcolor='blue')

    fig.tight_layout()
    fig.savefig(os.path.join(exp_path, name + '.png'))
    fig.clf()


def plot_adjacency(adjacency, gt_adjacency, exp_path, name=''):
    plt.figure(1)
    plt.subplot(121)
    plt.title("Model adj minus gt adj")
    plt.imshow(adjacency - gt_adjacency, vmin=-1, vmax=1)
    #for (j, i), mistake in np.ndenumerate(mistakes):
    #    if mistake:
    #        plt.text(i, j, ".", ha='center', va='center', color='r')

    plt.subplot(122)
    plt.title("Ground truth adjacency matrix")
    plt.imshow(gt_adjacency)

    plt.savefig(os.path.join(exp_path, 'adjacency' + name + '.png'))
    plt.clf()
    plt.close('all')


def plot_learning_curves(not_nlls, aug_lagrangians, aug_lagrangians_ma, aug_lagrangians_val, exp_path):
    fig, ax1 = plt.subplots()
    aug_lagrangians_val = np.array(aug_lagrangians_val)
    ax1.plot(range(len(not_nlls)), not_nlls, color="m", linewidth=1, label="Augmented Lagrangian minus NLL")
    ax1.plot(range(len(aug_lagrangians)), aug_lagrangians, color="k", linewidth=1, alpha=0.5)
    ax1.plot(range(len(aug_lagrangians_ma)), aug_lagrangians_ma, color="k", linewidth=1, label="Augmented Lagrangian")
    ax1.plot(aug_lagrangians_val[:,0], aug_lagrangians_val[:,1], color="r", linewidth=1,
             label="Augmented Lagrangian on validation set")

    ax1.set_xlabel("Iterations")
    ax1.set_ylim(bottom=0, top=2)

    ax1.legend()

    fig.tight_layout()
    fig.savefig(os.path.join(exp_path, 'learning-curves.png'))
    fig.clf()


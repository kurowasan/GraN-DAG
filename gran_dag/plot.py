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

import matplotlib

# To avoid displaying the figures
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns


def plot_weighted_adjacency(weighted_adjacency, gt_adjacency, exp_path, name="abs-weight-product", mus=None,
                            lambdas=None, iter=None, plotting_callback=None):
    """iter is useful to deal with jacobian, it will interpolate."""
    num_vars = weighted_adjacency.shape[1]
    max_value = 0
    fig, ax1 = plt.subplots()

    # Plot weight of incorrect edges
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

    # Plot weight of correct edges
    for i in range(num_vars):
        for j in range(num_vars):
            if gt_adjacency[i, j]:
                color = 'g'
            else:
                continue
            y = weighted_adjacency[:, i, j]
            num_iter = len(y)
            # plt.plot(range(len(weighted_adjacency[:, 0, 0])), y, color, alpha=0.1, linewidth=1)
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
    if plotting_callback is not None:
        plotting_callback("weighted_adjacency", fig)
    fig.savefig(os.path.join(exp_path, name + '.png'))
    fig.clf()


def plot_adjacency(adjacency, gt_adjacency, exp_path, name=''):
    plt.clf()
    f, (ax1, ax2, ax3) = plt.subplots(ncols=3, nrows=1)
    sns.heatmap(adjacency, ax=ax1, cbar=False, vmin=-1, vmax=1, cmap="Blues_r", xticklabels=False, yticklabels=False)
    sns.heatmap(gt_adjacency, ax=ax2, cbar=False, vmin=-1, vmax=1, cmap="Blues_r", xticklabels=False, yticklabels=False)
    sns.heatmap(adjacency - gt_adjacency, ax=ax3, cbar=False, vmin=-1, vmax=1, cmap="Blues_r", xticklabels=False,
                yticklabels=False)

    ax1.set_title("Learned")
    ax2.set_title("Ground truth")
    ax3.set_title("Learned - GT")

    ax1.set_aspect('equal', adjustable='box')
    ax2.set_aspect('equal', adjustable='box')
    ax3.set_aspect('equal', adjustable='box')

    plt.savefig(os.path.join(exp_path, 'adjacency' + name + '.png'))


def plot_learning_curves(not_nlls, aug_lagrangians, aug_lagrangians_ma, aug_lagrangians_val, nlls, nlls_val, exp_path):
    fig, ax1 = plt.subplots()
    aug_lagrangians_val = np.array(aug_lagrangians_val)
    ax1.plot(range(len(nlls)), nlls, color="orange", linewidth=1, label="NLL")
    ax1.plot(aug_lagrangians_val[:,0], nlls_val, color="blue", linewidth=1, label="NLL on validation set")
    print(len(nlls_val))
    ax1.plot(range(len(not_nlls)), not_nlls, color="m", linewidth=1, label="AL minus NLL")
    ax1.plot(range(len(aug_lagrangians)), aug_lagrangians, color="k", linewidth=1, alpha=0.5)
    ax1.plot(range(len(aug_lagrangians_ma)), aug_lagrangians_ma, color="k", linewidth=1, label="Augmented Lagrangian")
    ax1.plot(aug_lagrangians_val[:,0], aug_lagrangians_val[:,1], color="r", linewidth=1,
             label="Augmented Lagrangian on validation set")

    ax1.set_xlabel("Iterations")
    ax1.set_ylim(bottom=0, top=2)

    ax1.legend()

    fig.tight_layout()
    fig.savefig(os.path.join(exp_path, 'learning-curves.pdf'), bbox_inches="tight", padding=0)
    fig.clf()

def plot_learning_curves_retrain(losses, losses_val, nlls, nlls_val, exp_path):
    fig, ax1 = plt.subplots()
    losses_val = np.array(losses_val)
    ax1.plot(range(len(nlls)), nlls, color="orange", linewidth=1, label="NLL")
    ax1.plot(losses_val[:,0], nlls_val, color="blue", linewidth=1, label="NLL (val)")
    ax1.plot(range(len(losses)), losses, color="k", linewidth=1, label="NLL + REG")
    ax1.plot(losses_val[:,0], losses_val[:,1], color="r", linewidth=1, label="NLL + REG (val)")
    ax1.set_xlabel("Iterations")
    ax1.set_ylim(bottom=0, top=2)

    ax1.legend()

    fig.tight_layout()
    fig.savefig(os.path.join(exp_path, 'learning-curves.pdf'), bbox_inches="tight", padding=0)
    fig.clf()
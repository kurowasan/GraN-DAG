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
import sys
#sys.path.remove("/usr/local/lib/python3.7/dist-packages/cdt-0.5.2-py3.7.egg")
import cdt
import os
import time
import math
import copy

import numpy as np

np.set_printoptions(linewidth=200)
from sklearn.ensemble import ExtraTreesRegressor
from sklearn.feature_selection import SelectFromModel
import torch
import pandas as pd
from cdt.utils.R import RPackages, launch_R_script

from .dag_optim import compute_constraint, compute_jacobian_avg, is_acyclic
from .utils.metrics import edge_errors, shd as shd_metric
from .utils.penalty import compute_penalty, compute_group_lasso_l2_penalty
from .utils.save import dump, load, np_to_csv
from .utils.topo_sort import  generate_complete_dag
from .plot import plot_weighted_adjacency, plot_adjacency, plot_learning_curves, plot_learning_curves_retrain


EPSILON = 1e-8


def pns_(model_adj, train_data, test_data, num_neighbors, thresh):
    """Preliminary neighborhood selection"""
    x_train, _ = train_data.sample(train_data.num_samples)
    x_test, _ = test_data.sample(test_data.num_samples)
    x = np.concatenate([x_train.detach().cpu().numpy(), x_test.detach().cpu().numpy()], 0)

    num_samples = x.shape[0]
    num_nodes = x.shape[1]
    print("PNS: num samples = {}, num nodes = {}".format(num_samples, num_nodes))
    for node in range(num_nodes):
        print("PNS: node " + str(node))
        x_other = np.copy(x)
        x_other[:, node] = 0
        reg = ExtraTreesRegressor(n_estimators=500)
        reg = reg.fit(x_other, x[:, node])
        selected_reg = SelectFromModel(reg, threshold="{}*mean".format(thresh), prefit=True,
                                       max_features=num_neighbors)
        mask_selected = selected_reg.get_support(indices=False).astype(np.float)

        model_adj[:, node] *= mask_selected

    return model_adj


def pns(model, train_data, test_data, num_neighbors, thresh, exp_path, metrics_callback, plotting_callback):
    # Prepare path for saving results
    save_path = os.path.join(exp_path, "pns")
    if not os.path.exists(save_path):
        os.makedirs(save_path)

    # Check if already computed
    if os.path.exists(os.path.join(save_path, "DAG.npy")):
        print("pns already computed. Loading result from disk.")
        return load(save_path, "model.pkl")

    model_adj = model.adjacency.detach().cpu().numpy()
    time0 = time.time()
    model_adj = pns_(model_adj, train_data, test_data, num_neighbors, thresh)

    with torch.no_grad():
        model.adjacency.copy_(torch.Tensor(model_adj))

    timing = time.time() - time0
    print("PNS took {}s".format(timing))

    # save
    dump(model, save_path, 'model')
    dump(timing, save_path, 'timing')
    np.save(os.path.join(save_path, "DAG"), model.adjacency.detach().cpu().numpy())

    # plot
    plot_adjacency(model.adjacency.detach().cpu().numpy(), train_data.adjacency.detach().cpu().numpy(), save_path)

    return model


def train(model, gt_adjacency, train_data, test_data, opt, metrics_callback, plotting_callback):
    """
    Applying augmented Lagrangian to solve the continuous constrained problem.

    """

    # Prepare path for saving results
    save_path = os.path.join(opt.exp_path, "train")
    if not os.path.exists(save_path):
        os.makedirs(save_path)

    # Check if already computed
    if os.path.exists(os.path.join(save_path, "DAG.npy")):
        print("Train already computed. Loading result from disk.")
        return load(save_path, "model.pkl")

    time0 = time.time()

    # initialize stuff for learning loop
    aug_lagrangians = []
    aug_lagrangian_ma = [0.0] * (opt.num_train_iter + 1)
    aug_lagrangians_val = []
    grad_norms = []
    grad_norm_ma = [0.0] * (opt.num_train_iter + 1)
    if not opt.no_w_adjs_log:
        w_adjs = np.zeros((opt.num_train_iter, opt.num_vars, opt.num_vars), dtype=np.float32)
    hs = []
    not_nlls = []  # Augmented Lagrangrian minus (pseudo) NLL
    nlls = []  # NLL on train
    nlls_val = []  # NLL on validation
    delta_mu = np.inf

    # Augmented Lagrangian stuff
    mu = opt.mu_init
    lamb = opt.lambda_init
    mus = []
    lambdas = []

    if opt.optimizer == "sgd":
        optimizer = torch.optim.SGD(model.parameters(), lr=opt.lr)
    elif opt.optimizer == "rmsprop":
        optimizer = torch.optim.RMSprop(model.parameters(), lr=opt.lr)
    else:
        raise NotImplementedError("optimizer {} is not implemented".format(opt.optimizer))

    # Learning loop:
    for iter in range(opt.num_train_iter):
        # compute loss
        model.train()
        x, _ = train_data.sample(opt.train_batch_size)
        weights, biases, extra_params = model.get_parameters(mode="wbx")
        loss = - torch.mean(model.compute_log_likelihood(x, weights, biases, extra_params))
        nlls.append(loss.item())
        model.eval()

        # constraint related
        w_adj = model.get_w_adj()
        h = compute_constraint(model, w_adj)

        # compute augmented langrangian
        aug_lagrangian = loss + 0.5 * mu * h ** 2 + lamb * h

        # optimization step on augmented lagrangian
        optimizer.zero_grad()
        aug_lagrangian.backward()
        optimizer.step()

        # clamp edges
        if opt.edge_clamp_range != 0:
            with torch.no_grad():
                to_keep = (w_adj > opt.edge_clamp_range).type(torch.Tensor)
                model.adjacency *= to_keep

        # logging
        if not opt.no_w_adjs_log:
            w_adjs[iter, :, :] = w_adj.detach().cpu().numpy().astype(np.float32)
        mus.append(mu)
        lambdas.append(lamb)
        not_nlls.append(0.5 * mu * h.item() ** 2 + lamb * h.item())

        # compute augmented lagrangian moving average
        aug_lagrangians.append(aug_lagrangian.item())
        aug_lagrangian_ma[iter + 1] = aug_lagrangian_ma[iter] + 0.01 * (aug_lagrangian.item() - aug_lagrangian_ma[iter])
        grad_norms.append(model.get_grad_norm("wbx").item())
        grad_norm_ma[iter + 1] = grad_norm_ma[iter] + 0.01 * (grad_norms[-1] - grad_norm_ma[iter])

        # compute loss on whole validation set
        if iter % opt.stop_crit_win == 0:
            with torch.no_grad():
                x, _ = test_data.sample(test_data.num_samples)
                loss_val = - torch.mean(model.compute_log_likelihood(x, weights, biases, extra_params)).item()
                nlls_val.append(loss_val)
                aug_lagrangians_val.append([iter, loss_val + not_nlls[-1]])

        # compute delta for lambda
        if iter >= 2 * opt.stop_crit_win and iter % (2 * opt.stop_crit_win) == 0:
            t0, t_half, t1 = aug_lagrangians_val[-3][1], aug_lagrangians_val[-2][1], aug_lagrangians_val[-1][1]

            # if the validation loss went up and down, do not update lagrangian and penalty coefficients.
            if not (min(t0, t1) < t_half < max(t0, t1)):
                delta_lambda = -np.inf
            else:
                delta_lambda = (t1 - t0) / opt.stop_crit_win
        else:
            delta_lambda = -np.inf  # do not update lambda nor mu

        # log metrics
        if iter % 100 == 0:
            print("Iteration:", iter)
            if opt.num_vars <= 5:
                print("    w_adj:\n", w_adj.detach().cpu().numpy())
                print("    current_adjacency:\n", model.adjacency.detach().cpu().numpy())
                print("    gt_adjacency:\n", gt_adjacency)

            # XXX: just for debug (graph metrics in real-time)
            with torch.no_grad():
                to_keep = (model.get_w_adj() > opt.edge_clamp_range).type(torch.Tensor)
                current_adj = model.adjacency * to_keep
                current_adj = current_adj.cpu().numpy()
                fn, fp, rev = edge_errors(current_adj, gt_adjacency)

            metrics_callback(stage="train", step=iter,
                             metrics={"aug-lagrangian": aug_lagrangian.item(),
                                      "aug-lagrangian-moving-avg": aug_lagrangian_ma[iter + 1],
                                      "aug-lagrangian-val": aug_lagrangians_val[-1][1],
                                      "nll": nlls[-1],
                                      "nll-val": nlls_val[-1],
                                      "nll-gap": nlls_val[-1] - nlls[-1],
                                      "grad-norm-moving-average": grad_norm_ma[iter + 1],
                                      "delta_lambda": delta_lambda,
                                      "omega_lambda": opt.omega_lambda,
                                      "delta_mu": delta_mu,
                                      "omega_mu": opt.omega_mu,
                                      "constraint_violation": h.item(),
                                      "mu": mu,
                                      "lambda": lamb,
                                      "initial_lr": opt.lr,
                                      "current_lr": opt.lr,
                                      "fn_edges": fn,
                                      "fp_edges": fp,
                                      "rev_edges": rev,
                                      "shd": fn + fp + rev
                                      })

        # plot
        if iter % opt.plot_freq == 0:
            if not opt.no_w_adjs_log:
                plot_weighted_adjacency(w_adjs[:iter + 1], gt_adjacency, opt.exp_path,
                                        name="w_adj", mus=mus, lambdas=lambdas,
                                        plotting_callback=plotting_callback)
            plot_adjacency(model.adjacency.detach().cpu().numpy(), gt_adjacency, opt.exp_path)
            plot_learning_curves(not_nlls, aug_lagrangians, aug_lagrangian_ma[:iter], aug_lagrangians_val, nlls,
                                 nlls_val, opt.exp_path)

        # Does the augmented lagrangian converged?
        if h > opt.h_threshold:
            # if we have found a stationary point of the augmented loss
            if abs(delta_lambda) < opt.omega_lambda or delta_lambda > 0:
                lamb += mu * h.item()
                print("Updated lambda to {}".format(lamb))

                # Did the constraint improve sufficiently?
                hs.append(h.item())
                if len(hs) >= 2:
                    if hs[-1] > hs[-2] * opt.omega_mu:
                        mu *= 10
                        print("Updated mu to {}".format(mu))

                # little hack to make sure the moving average is going down.
                with torch.no_grad():
                    gap_in_not_nll = 0.5 * mu * h.item() ** 2 + lamb * h.item() - not_nlls[-1]
                    aug_lagrangian_ma[iter + 1] += gap_in_not_nll
                    aug_lagrangians_val[-1][1] += gap_in_not_nll

                if opt.optimizer == "rmsprop":
                    optimizer = torch.optim.RMSprop(model.parameters(), lr=opt.lr_reinit)
                else:
                    optimizer = torch.optim.SGD(model.parameters(), lr=opt.lr_reinit)

        else:
            timing = time.time() - time0

            # Final clamping of all edges == 0
            with torch.no_grad():
                to_keep = (w_adj > 0).type(torch.Tensor)
                model.adjacency *= to_keep

            # compute nll on train and validation set
            weights, biases, extra_params = model.get_parameters(mode="wbx")
            x, _ = train_data.sample(train_data.num_samples)
            # Since we do not have a DAG yet, this is not really a negative log likelihood.
            nll_train = - torch.mean(model.compute_log_likelihood(x, weights, biases, extra_params))

            x, _ = test_data.sample(test_data.num_samples)
            nll_val = - torch.mean(model.compute_log_likelihood(x, weights, biases, extra_params))

            # Save
            if not opt.no_w_adjs_log:
                w_adjs = w_adjs[:iter]
            dump(model, save_path, 'model')
            dump(opt, save_path, 'opt')
            if opt.num_vars <= 50 and not opt.no_w_adjs_log:
                dump(w_adjs, save_path, 'w_adjs')
            dump(nll_train, save_path, 'pseudo-nll-train')
            dump(nll_val, save_path, 'pseudo-nll-val')
            dump(not_nlls, save_path, 'not-nlls')
            dump(aug_lagrangians, save_path, 'aug-lagrangians')
            dump(aug_lagrangian_ma[:iter], save_path, 'aug-lagrangian-ma')
            dump(aug_lagrangians_val, save_path, 'aug-lagrangians-val')
            dump(grad_norms, save_path, 'grad-norms')
            dump(grad_norm_ma[:iter], save_path, 'grad-norm-ma')
            dump(timing, save_path, 'timing')
            np.save(os.path.join(save_path, "DAG"), model.adjacency.detach().cpu().numpy())

            # plot
            if not opt.no_w_adjs_log:
                plot_weighted_adjacency(w_adjs, gt_adjacency, save_path,
                                        name="w_adj", mus=mus, lambdas=lambdas)
            plot_adjacency(model.adjacency.detach().cpu().numpy(), gt_adjacency, save_path)
            plot_learning_curves(not_nlls, aug_lagrangians, aug_lagrangian_ma[:iter], aug_lagrangians_val, nlls,
                                 nlls_val, save_path)

            return model


def to_dag(model, train_data, test_data, opt, metrics_callback, plotting_callback, stage_name="to-dag"):
    """
    1- If some entries of A_\phi == 0, also mask them (This can happen with stochastic proximal gradient descent)
    2- Remove edges (from weaker to stronger) until a DAG is obtained.

    """
    # Prepare path for saving results
    save_path = os.path.join(opt.exp_path, stage_name)
    if not os.path.exists(save_path):
        os.makedirs(save_path)

    # Check if already computed
    if os.path.exists(os.path.join(save_path, "DAG.npy")):
        print(stage_name + " already computed. Loading result from disk.")
        return load(save_path, "model.pkl")

    model.eval()
    time0 = time.time()

    if opt.jac_thresh:
        A = compute_jacobian_avg(model, train_data, train_data.num_samples).t()
    else:
        A = model.get_w_adj()
    A = A.detach().cpu().numpy()

    with torch.no_grad():
        # Find the smallest threshold that removes all cycle-inducing edges
        thresholds = np.unique(A)
        for step, t in enumerate(thresholds):
            print("Edges/thresh", model.adjacency.sum(), t)
            to_keep = torch.Tensor(A > t + EPSILON)
            new_adj = model.adjacency * to_keep

            metrics_callback(stage=stage_name, step=step,
                             metrics={"threshold": t, "edges": new_adj.sum().item()},
                             throttle=False)

            if is_acyclic(new_adj):
                model.adjacency.copy_(new_adj)
                break

    # evaluate on validation set
    x, _ = test_data.sample(test_data.num_samples)
    weights, biases, extra_params = model.get_parameters(mode="wbx")
    nll_val = - torch.mean(model.compute_log_likelihood(x, weights, biases, extra_params)).item()

    # Compute SHD and SID metrics
    pred_adj_ = model.adjacency.detach().cpu().numpy()
    train_adj_ = train_data.adjacency.detach().cpu().numpy()
    sid = float(cdt.metrics.SID(target=train_adj_, pred=pred_adj_))
    shd = float(shd_metric(pred_adj_, train_adj_))
    shd_cpdag = float(cdt.metrics.SHD_CPDAG(target=train_adj_, pred=pred_adj_))
    fn, fp, rev = edge_errors(pred_adj_, train_adj_)
    del train_adj_, pred_adj_

    metrics_callback(stage=stage_name, step=0,
                     metrics={"nll_val": nll_val, "sid": sid, "shd": shd, "shd_cpdag": shd_cpdag, "fn": fn, "fp": fp,
                              "rev": rev},
                     throttle=False)

    timing = time.time() - time0

    # Save
    dump(model, save_path, 'model')
    dump(opt, save_path, 'opt')
    dump(timing, save_path, 'timing')
    dump(nll_val, save_path, "nll-val", txt=True)
    np.save(os.path.join(save_path, "DAG"), model.adjacency.detach().cpu().numpy())

    # plot adjacency
    plot_adjacency(model.adjacency.detach().cpu().numpy(), train_data.adjacency.detach().cpu().numpy(), save_path)

    return model


def cam_pruning_(model_adj, train_data, test_data, cutoff, save_path, verbose=False):
    # convert numpy data to csv, so R can access them
    data_np = np.concatenate([train_data.dataset.detach().cpu().numpy(), test_data.dataset.detach().cpu().numpy()], 0)
    data_csv_path = np_to_csv(data_np, save_path)
    dag_csv_path = np_to_csv(model_adj, save_path)

    #dag_pruned = cam_pruning(path_data, path_dag, cutoff, verbose)
    if not RPackages.CAM:
        raise ImportError("R Package CAM is not available.")

    arguments = dict()
    arguments['{PATH_DATA}'] = data_csv_path
    arguments['{PATH_DAG}'] = dag_csv_path
    arguments['{PATH_RESULTS}'] = os.path.join(save_path, "results.csv")
    arguments['{CUTOFF}'] = str(cutoff)

    if verbose:
        arguments['{VERBOSE}'] = "TRUE"
    else:
        arguments['{VERBOSE}'] = "FALSE"

    def retrieve_result():
        return pd.read_csv(arguments['{PATH_RESULTS}']).values

    dag_pruned = launch_R_script("{}/utils/cam_pruning.R".format(os.path.dirname(os.path.realpath(__file__))),
                                     arguments, output_function=retrieve_result)

    # remove the temporary csv files
    os.remove(data_csv_path)
    os.remove(dag_csv_path)

    return dag_pruned


def cam_pruning(model, train_data, test_data, opt, metrics_callback, plotting_callback, cutoff=0.001, verbose=False):
    """Execute CAM pruning on a given model and datasets"""
    time0 = time.time()
    model.eval()

    # Prepare path for saving results
    stage_name = "cam-pruning/cutoff_%.6f" % cutoff
    save_path = os.path.join(opt.exp_path, stage_name)
    if not os.path.exists(save_path):
        os.makedirs(save_path)

    # Check if already computed
    if os.path.exists(os.path.join(save_path, "DAG.npy")):
        print(stage_name, "already computed. Loading result from disk.")
        return load(save_path, "model.pkl")

    model_adj = model.adjacency.detach().cpu().numpy()

    dag_pruned = cam_pruning_(model_adj, train_data, test_data, cutoff, save_path, verbose)

    # set new adjacency matrix to model
    model.adjacency.copy_(torch.as_tensor(dag_pruned).type(torch.Tensor))

    # evaluate on validation set
    x, _ = test_data.sample(test_data.num_samples)
    weights, biases, extra_params = model.get_parameters(mode="wbx")
    nll_val = - torch.mean(model.compute_log_likelihood(x, weights, biases, extra_params)).item()

    # Compute SHD and SID metrics
    pred_adj_ = model.adjacency.detach().cpu().numpy()
    train_adj_ = train_data.adjacency.detach().cpu().numpy()
    sid = float(cdt.metrics.SID(target=train_adj_, pred=pred_adj_))
    shd = float(shd_metric(pred_adj_, train_adj_))
    shd_cpdag = float(cdt.metrics.SHD_CPDAG(target=train_adj_, pred=pred_adj_))
    fn, fp, rev = edge_errors(pred_adj_, train_adj_)
    del train_adj_, pred_adj_

    metrics_callback(stage=stage_name, step=0,
                     metrics={"nll_val": nll_val, "sid": sid, "shd": shd, "shd_cpdag": shd_cpdag, "fn": fn, "fp": fp,
                              "rev": rev}, throttle=False)

    timing = time.time() - time0

    # save
    dump(model, save_path, 'model')
    dump(opt, save_path, 'opt')
    dump(timing, save_path, 'timing')
    dump(nll_val, save_path, "nll-val", txt=True)
    np.save(os.path.join(save_path, "DAG"), model.adjacency.detach().cpu().numpy())

    # plot adjacency
    plot_adjacency(model.adjacency.detach().cpu().numpy(), train_data.adjacency.detach().cpu().numpy(), save_path)

    return model


def retrain(model, train_data, test_data, dag_folder, opt, metrics_callback, plotting_callback):
    """
    Retrain a model which is already DAG

    """
    # Prepare path for saving results
    stage_name = "retrain_{}".format(dag_folder)
    save_path = os.path.join(opt.exp_path, stage_name)
    if not os.path.exists(save_path):
        os.makedirs(save_path)

    # Check if already computed
    if os.path.exists(os.path.join(save_path, "best-model.pkl")):
        print(stage_name, "already computed. Loading result from disk.")
        return load(save_path, "best-model.pkl")

    time0 = time.time()

    # initialize stuff for learning loop
    nlls = []
    nlls_val = []
    losses = []
    losses_val = []
    grad_norms = []
    grad_norm_ma = [0.0] * (opt.num_train_iter + 1)

    # early stopping stuff
    best_model = copy.deepcopy(model)
    best_nll_val = np.inf
    patience = opt.patience

    if opt.optimizer == "sgd":
        optimizer = torch.optim.SGD(model.parameters(), lr=opt.lr)
    elif opt.optimizer == "rmsprop":
        optimizer = torch.optim.RMSprop(model.parameters(), lr=opt.lr)
    else:
        raise NotImplementedError("optimizer {} is not implemented".format(opt.optimizer))

    # Learning loop:
    for iter in range(opt.num_train_iter):
        # compute loss
        model.train()
        x, _ = train_data.sample(opt.train_batch_size)
        weights, biases, extra_params = model.get_parameters(mode="wbx")
        nll = - torch.mean(model.compute_log_likelihood(x, weights, biases, extra_params))
        nlls.append(nll.item())
        model.eval()

        loss = nll

        # optimization step on augmented lagrangian
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()

        # compute augmented lagrangian moving average
        losses.append(loss.item())
        grad_norms.append(model.get_grad_norm("wbx").item())
        grad_norm_ma[iter + 1] = grad_norm_ma[iter] + 0.01 * (grad_norms[-1] - grad_norm_ma[iter])

        # compute loss on whole validation set
        if iter % 1000 == 0:
            with torch.no_grad():
                # import ipdb; ipdb.set_trace(context=5)
                x, _ = test_data.sample(test_data.num_samples)
                nll_val = - torch.mean(model.compute_log_likelihood(x, weights, biases, extra_params)).item()
                nlls_val.append(nll_val)
                losses_val.append([iter, nll_val])

                # nll_val the best?
                if nll_val < best_nll_val:
                    best_nll_val = nll_val
                    patience = opt.patience
                    best_model = copy.deepcopy(model)
                else:
                    patience -= 1

        # log metrics
        if iter % 100 == 0:
            print("Iteration:", iter)
            metrics_callback(stage=stage_name, step=iter,
                             metrics={"loss": loss.item(),
                                      "loss-val": losses_val[-1][1],
                                      "nll": nlls[-1],
                                      "nll-val": nlls_val[-1],
                                      "grad-norm-moving-average": grad_norm_ma[iter + 1],
                                      "patience": patience,
                                      "best-nll-val": best_nll_val})

        # plot
        if iter % opt.plot_freq == 0:
            plot_learning_curves_retrain(losses, losses_val, nlls, nlls_val, save_path)

        # Have we converged?
        if patience == 0:
            timing = time.time() - time0

            # save
            # call to-dag in order to mask the edges i,j for which (A_\phi)_ij == 0
            best_model = to_dag(best_model, train_data, test_data, opt, metrics_callback, plotting_callback,
                                stage_name="final_to-dag")
            dump(best_model, save_path, 'best-model')
            dump(best_nll_val, save_path, 'best-nll-val', txt=True)
            dump(opt, save_path, 'opt')
            dump(nlls, save_path, 'nlls-train')
            dump(nlls_val, save_path, 'nlls-val')
            dump(losses, save_path, 'losses')
            dump(losses_val, save_path, 'losses-val')
            dump(grad_norms, save_path, 'grad-norms')
            dump(grad_norm_ma[:iter], save_path, 'grad-norm-ma')
            dump(timing, save_path, 'timing')

            # plot
            plot_learning_curves_retrain(losses, losses_val, nlls, nlls_val, save_path)

            return model

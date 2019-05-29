import os
import time

import numpy as np
np.set_printoptions(linewidth=200)
from sklearn.ensemble import ExtraTreesRegressor
from sklearn.feature_selection import SelectFromModel
import torch

from dag_optim import compute_constraint, compute_prod, compute_jacobian_avg, is_acyclic
from utils import load, dump
from plot import plot_weighted_adjacency, plot_adjacency, plot_learning_curves


def pns(model, train_data, test_data, num_neighbors, thresh, exp_path):
    """Preliminary neighborhood selection"""
    time0 = time.time()
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
        selected_reg = SelectFromModel(reg, threshold="{}*mean".format(thresh) , prefit=True,
                                       max_features=num_neighbors)
        mask_selected = torch.Tensor(selected_reg.get_support(indices=False).astype(np.float))

        with torch.no_grad():
            model.adjacency[:, node] *= mask_selected

    timing = time.time() - time0
    print("PNS took {}s".format(timing))
    # save
    save_path = os.path.join(exp_path, "pns")
    if not os.path.exists(save_path):
        os.makedirs(save_path)
    dump(model, save_path, 'model')
    dump(timing, save_path, 'timing')
    np.save(os.path.join(save_path, "DAG"), model.adjacency.detach().cpu().numpy())

    # plot
    plot_adjacency(model.adjacency.detach().cpu().numpy(), train_data.adjacency.detach().cpu().numpy(), save_path)

    return model

def train(model, gt_adjacency, train_data, test_data, opt):
    """Applying augmented Lagrangian to solve the continuous constrained problem."""
    time0 = time.time()
    # Initialize stuff for learning loop
    optimizer = torch.optim.RMSprop(model.parameters(), lr=opt.lr)
    aug_lagrangians = []
    aug_lagrangian_ma = [0.0] * (opt.num_train_iter + 1)
    aug_lagrangians_val = []
    grad_norms = []
    grad_norm_ma = [0.0] * (opt.num_train_iter + 1)
    prods = np.zeros((opt.num_train_iter, opt.num_vars, opt.num_vars), dtype=np.float32)
    h = np.inf  # constraint violation function
    hs = []
    not_nlls = []  # Augmented Lagrangrian minus (pseudo) NLL
    delta_mu = np.inf

    # Augmented Lagrangian stuff
    mu = opt.mu_init
    lamb = opt.lambda_init

    mus = []
    lambdas = []

    # Learning loop:
    for iter in range(opt.num_train_iter):
        # compute loss
        x, _ = train_data.sample(opt.train_batch_size)
        weights, biases, extra_params = model.get_parameters(mode="wbx")
        loss = - torch.mean(model.compute_log_likelihood(x, weights, biases, extra_params))

        # constraint related
        prod = compute_prod(model, norm=opt.norm_prod)
        if h != 0:
            h = compute_constraint(model, prod)

        # compute augmented langrangian
        aug_lagrangian = loss + 0.5 * mu * h ** 2 + lamb * h

        # optimization step on augmented lagrangian
        optimizer.zero_grad()
        aug_lagrangian.backward()
        optimizer.step()

        # clamp edges
        if opt.edge_clamp_range is not None:
            with torch.no_grad():
                to_keep = (prod > opt.edge_clamp_range).t().type(torch.Tensor)
                model.adjacency *= to_keep

        # logging
        prods[iter, :, :] = prod.detach().cpu().numpy().astype(np.float32)
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
                #import ipdb; ipdb.set_trace(context=5)
                x, _ = test_data.sample(test_data.num_samples)
                loss_val = - torch.mean(model.compute_log_likelihood(x, weights, biases, extra_params)).item()
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
            delta_lambda = -np.inf # do not update lambda nor mu

        # log metrics
        if iter % 100 == 0:
            print("Iteration:", iter)
            if opt.num_vars <= 5:
                print("    prod.T:\n", prod.t().detach().cpu().numpy())
                print("    current_adjacency:\n", model.adjacency.detach().cpu().numpy())
                print("    gt_adjacency:\n", gt_adjacency)
            print("    aug-lagrangian:", aug_lagrangian.item())
            print("    aug-lagrangian-moving-avg:", aug_lagrangian_ma[iter + 1])
            print("    aug-lagrangian-val:", aug_lagrangians_val[-1][1])
            print("    grad-norm-moving-average:", grad_norm_ma[iter + 1])
            print("    delta_lambda:", delta_lambda)
            print("    omega_lambda:", opt.omega_lambda)
            print("    delta_mu:", delta_mu)
            print("    omega_mu:", opt.omega_mu)
            print("    constraint violation:", h.item())
            print("    mu:", mu)
            print("    lambda:", lamb)

        # plot
        if iter % opt.plot_freq == 0:
            plot_weighted_adjacency(prods[:iter + 1].swapaxes(1, 2), gt_adjacency, opt.exp_path,
                                    name="abs-weight-product", mus=mus, lambdas=lambdas)
            plot_adjacency(model.adjacency.detach().cpu().numpy(), gt_adjacency, opt.exp_path)
            plot_learning_curves(not_nlls, aug_lagrangians, aug_lagrangian_ma[:iter], aug_lagrangians_val, opt.exp_path)

        # Does the augmented lagrangian converged?
        if h > opt.h_threshold:
            # if we have found a stationary point of the augmented loss
            if abs(delta_lambda) < opt.omega_lambda or delta_lambda > 0:
                lamb += mu * h.item()

                # Did the constraint improve sufficiently?
                hs.append(h.item())
                if len(hs) >= 2:
                    if hs[-1] > hs[-2] * opt.omega_mu:
                        mu *= 10

                # little hack to make sure the moving average is going down.
                with torch.no_grad():
                    gap_in_not_nll = 0.5 * mu * h.item() ** 2 + lamb * h.item() - not_nlls[-1]
                    aug_lagrangian_ma[iter + 1] += gap_in_not_nll
                    aug_lagrangians_val[-1][1] += gap_in_not_nll

                # reinitializing the optimizer
                optimizer = torch.optim.RMSprop(model.parameters(), lr=opt.lr_reinit)

        else:
            timing = time.time() - time0

            # compute nll on train and validation set
            weights, biases, extra_params = model.get_parameters(mode="wbx")
            x, _ = train_data.sample(train_data.num_samples)
            nll_train = - torch.mean(model.compute_log_likelihood(x, weights, biases, extra_params))

            x, _ = test_data.sample(test_data.num_samples)
            nll_val = - torch.mean(model.compute_log_likelihood(x, weights, biases, extra_params))

            # save
            prods = prods[:iter]
            save_path = os.path.join(opt.exp_path, "train")
            if not os.path.exists(save_path):
                os.makedirs(save_path)
            dump(model, save_path, 'model')
            dump(opt, save_path, 'opt')
            if opt.num_vars <=50:
                dump(prods, save_path, 'prods')
            dump(nll_train, save_path, 'nll-train')
            dump(nll_val, save_path, 'nll-val')
            dump(not_nlls, save_path, 'not-nlls')
            dump(aug_lagrangians, save_path, 'aug-lagrangians')
            dump(aug_lagrangian_ma[:iter], save_path, 'aug-lagrangian-ma')
            dump(aug_lagrangians_val, save_path, 'aug-lagrangians-val')
            dump(grad_norms, save_path, 'grad-norms')
            dump(grad_norm_ma[:iter], save_path, 'grad-norm-ma')
            dump(timing, save_path, 'timing')
            np.save(os.path.join(save_path, "DAG"), model.adjacency.detach().cpu().numpy())

            # plot
            plot_weighted_adjacency(prods.swapaxes(1, 2), gt_adjacency, save_path, mus=mus, lambdas=lambdas)
            plot_adjacency(model.adjacency.detach().cpu().numpy(), gt_adjacency, save_path)
            plot_learning_curves(not_nlls, aug_lagrangians, aug_lagrangian_ma[:iter], aug_lagrangians_val, save_path)

            return model


def to_dag(model, train_data, test_data, opt):
    """Threshold the edges (using the mask model.adjacency) to make sure a DAG is obtained"""
    time0 = time.time()
    # constraint is not fully satisfied yet, we remove edges until we have a DAG
    if opt.jac_thresh:
        avg_jac = compute_jacobian_avg(model, train_data, train_data.num_samples)
        order = torch.argsort(avg_jac.t().contiguous().view(-1))
    else:
        prod = compute_prod(model, norm=opt.norm_prod)
        order = torch.argsort(prod.t().contiguous().view(-1))
    with torch.no_grad():
        i = 0
        while not is_acyclic(model.adjacency.detach().cpu().numpy()):
            model.adjacency.reshape(-1)[order[i]] = 0
            i += 1

    # evaluate on validation set
    x, _ = test_data.sample(test_data.num_samples)
    weights, biases, extra_params = model.get_parameters(mode="wbx")
    nll_val = - torch.mean(model.compute_log_likelihood(x, weights, biases, extra_params)).item()

    timing = time.time() - time0

    # save
    save_path = os.path.join(opt.exp_path, "to-dag")
    if not os.path.exists(save_path):
        os.makedirs(save_path)
    dump(model, save_path, 'model')
    dump(opt, save_path, 'opt')
    dump(timing, save_path, 'timing')
    dump(nll_val, save_path, "best-val", txt=True)
    np.save(os.path.join(save_path, "DAG"), model.adjacency.detach().cpu().numpy())

    # plot adjacency
    plot_adjacency(model.adjacency.detach().cpu().numpy(), train_data.adjacency.detach().cpu().numpy(), save_path)

    return model


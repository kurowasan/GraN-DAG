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
import timeit

import numpy as np
from scipy.linalg import expm
from scipy.special import comb
import torch


class TrExpScipy(torch.autograd.Function):
    """
    autograd.Function to compute trace of an exponential of a matrix
    """

    @staticmethod
    def forward(ctx, input):
        with torch.no_grad():
            # send tensor to cpu in numpy format and compute expm using scipy
            expm_input = expm(input.detach().cpu().numpy())
            # transform back into a tensor
            expm_input = torch.as_tensor(expm_input)
            if input.is_cuda:
                # expm_input = expm_input.cuda()
                assert expm_input.is_cuda
            # save expm_input to use in backward
            ctx.save_for_backward(expm_input)

            # return the trace
            return torch.trace(expm_input)

    @staticmethod
    def backward(ctx, grad_output):
        with torch.no_grad():
            expm_input, = ctx.saved_tensors
            return expm_input.t() * grad_output


def compute_constraint(model, w_adj):
    assert (w_adj >= 0).detach().cpu().numpy().all()
    h = TrExpScipy.apply(w_adj) - model.num_vars
    return h


def is_acyclic(adjacency):
    prod = np.eye(adjacency.shape[0])
    for _ in range(1, adjacency.shape[0] + 1):
        prod = np.matmul(adjacency, prod)
        if np.trace(prod) != 0: return False
    return True


def compute_A_phi(model, norm="none", square=False):
    weights = model.get_parameters(mode='w')[0]
    prod = torch.eye(model.num_vars)
    if norm != "none":
        prod_norm = torch.eye(model.num_vars)
    for i, w in enumerate(weights):
        if square:
            w = w ** 2
        else:
            w = torch.abs(w)
        if i == 0:
            prod = torch.einsum("tij,ljt,jk->tik", w, model.adjacency.unsqueeze(0), prod)
            if norm != "none":
                tmp = 1. - torch.eye(model.num_vars).unsqueeze(0)
                prod_norm = torch.einsum("tij,ljt,jk->tik", torch.ones_like(w).detach(), tmp, prod_norm)
        else:
            prod = torch.einsum("tij,tjk->tik", w, prod)
            if norm != "none":
                prod_norm = torch.einsum("tij,tjk->tik", torch.ones_like(w).detach(), prod_norm)

    # sum over density parameter axis
    prod = torch.sum(prod, 1)
    if norm == "paths":
        prod_norm = torch.sum(prod_norm, 1)
        denominator = prod_norm + torch.eye(model.num_vars)  # avoid / 0 on diagonal
        return (prod / denominator).t()
    elif norm == "none":
        return prod.t()
    else:
        raise NotImplementedError


def compute_jacobian_avg(model, data_manager, batch_size):
    jac_avg = torch.zeros(model.num_vars, model.num_vars)

    # sample
    x, do_mask = data_manager.sample(batch_size)
    x.requires_grad = True

    # compute loss
    weights, biases, extra_params = model.get_parameters(mode="wbx")
    log_probs = model.compute_log_likelihood(x, weights, biases, extra_params, detach=True)
    log_probs = torch.unbind(log_probs, 1)

    # compute jacobian of the loss
    for i in range(model.num_vars):
        tmp = torch.autograd.grad(log_probs[i], x, retain_graph=True, grad_outputs=torch.ones(batch_size))[0]
        jac_avg[i, :] = torch.abs(tmp).mean(0)

    return jac_avg


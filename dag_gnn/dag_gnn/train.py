
''''
Main function for traininng dag_gnn

'''


from __future__ import division
from __future__ import print_function

import time
import argparse
import pickle
import os
import datetime
import copy

# import torch
import torch.optim as optim
from torch.optim import lr_scheduler
import math

import numpy as np
from .utils import *
from .modules import *

# compute constraint h(A) value
def _h_A(A, m):
    expm_A = matrix_poly(A*A, m)
    h_A = torch.trace(expm_A) - m
    return h_A

prox_plus = torch.nn.Threshold(0.,0.)

def stau(w, tau):
    w1 = prox_plus(torch.abs(w)-tau)
    return torch.sign(w)*w1


def update_optimizer(optimizer, original_lr, c_A):
    '''related LR to c_A, whenever c_A gets big, reduce LR proportionally'''
    MAX_LR = 1e-2
    MIN_LR = 1e-4

    estimated_lr = original_lr / (math.log10(c_A) + 1e-10)
    if estimated_lr > MAX_LR:
        lr = MAX_LR
    elif estimated_lr < MIN_LR:
        lr = MIN_LR
    else:
        lr = estimated_lr

    # set LR
    for parame_group in optimizer.param_groups:
        parame_group['lr'] = lr

    return optimizer, lr


def evaluate_score(encoder, decoder, valid_loader, w_threshold, rel_rec, rel_send, num_nodes, args):
    encoder.eval()
    decoder.eval()

    losses = []
    with torch.no_grad():
        for batch_idx, (data, relations) in enumerate(valid_loader):

            if args.cuda:
                data, relations = data.cuda(), relations.cuda()

            # reshape data
            relations = relations.unsqueeze(2)

            enc_x, logits, origin_A, adj_A_tilt_encoder, z_gap, z_positive, myA, Wa = encoder(data, rel_rec, rel_send)  # logits is of size: [num_sims, z_dims]
            edges = logits

            # Thresholding
            graph = torch.where(origin_A < w_threshold, torch.tensor(0.), origin_A)

            dec_x, output, adj_A_tilt_decoder = decoder(data, edges, num_nodes * args.x_dims, rel_rec, rel_send, graph, adj_A_tilt_encoder, Wa)

            if torch.sum(output != output):
                print('nan error\n')

            target = data
            preds = output
            variance = 0.

            # reconstruction accuracy loss
            loss_nll = nll_gaussian(preds, target, variance)

            # KL loss
            loss_kl = kl_gaussian_sem(logits)

            # ELBO loss:
            loss = loss_kl + loss_nll
            losses.append(loss)

    return -np.mean(losses)


def dag_gnn(args, train_data, gt_adjacency, mask_A):
    args.cuda = not args.no_cuda and torch.cuda.is_available()
    args.factor = not args.no_factor
    args.save_folder = args.exp_path

    if args.dynamic_graph:
        print("Testing with dynamically re-computed graph.")

    # Save model and meta-data. Always saves in a new sub-folder.
    if args.save_folder:
        exp_counter = 0
        now = datetime.datetime.now()
        timestamp = now.isoformat()
        save_folder = '{}/exp{}/'.format(args.save_folder, timestamp)
        # safe_name = save_folder.text.replace('/', '_')
        os.makedirs(save_folder)
        meta_file = os.path.join(save_folder, 'metadata.pkl')
        encoder_file = os.path.join(save_folder, 'encoder.pt')
        decoder_file = os.path.join(save_folder, 'decoder.pt')

        log_file = os.path.join(save_folder, 'log.txt')
        log = open(log_file, 'w')

        pickle.dump({'args': args}, open(meta_file, "wb"))
    else:
        print("WARNING: No save_folder provided!" +
              "Testing (within this script) will throw an error.")

    # ================================================
    # get data: experiments = {synthetic SEM, ALARM}
    # ================================================
    train_loader, ground_truth_G = load_numpy_data(args, train_data, gt_adjacency)

    # ===================================
    # load modules
    # ===================================
    # Generate off-diagonal interaction graph
    off_diag = np.ones([args.data_variable_size, args.data_variable_size]) - np.eye(args.data_variable_size)

    rel_rec = np.array(encode_onehot(np.where(off_diag)[1]), dtype=np.float64)
    rel_send = np.array(encode_onehot(np.where(off_diag)[0]), dtype=np.float64)
    rel_rec = torch.DoubleTensor(rel_rec)
    rel_send = torch.DoubleTensor(rel_send)

    # add adjacency matrix A
    num_nodes = args.data_variable_size
    adj_A = np.zeros((num_nodes, num_nodes))

    if args.encoder == 'mlp':
        encoder = MLPEncoder(args.data_variable_size * args.x_dims, args.x_dims, args.encoder_hidden,
                             int(args.z_dims), adj_A, mask_A,
                             batch_size=args.batch_size,
                             do_prob=args.encoder_dropout, factor=args.factor).double()
    elif args.encoder == 'sem':
        encoder = SEMEncoder(args.data_variable_size * args.x_dims, args.encoder_hidden,
                             int(args.z_dims), adj_A,
                             batch_size=args.batch_size,
                             do_prob=args.encoder_dropout, factor=args.factor).double()

    if args.decoder == 'mlp':
        decoder = MLPDecoder(args.data_variable_size * args.x_dims,
                             args.z_dims, args.x_dims, encoder,
                             data_variable_size=args.data_variable_size,
                             batch_size=args.batch_size,
                             n_hid=args.decoder_hidden,
                             do_prob=args.decoder_dropout).double()
    elif args.decoder == 'sem':
        decoder = SEMDecoder(args.data_variable_size * args.x_dims,
                             args.z_dims, 2, encoder,
                             data_variable_size=args.data_variable_size,
                             batch_size=args.batch_size,
                             n_hid=args.decoder_hidden,
                             do_prob=args.decoder_dropout).double()

    if args.load_folder:
        encoder_file = os.path.join(args.load_folder, 'encoder.pt')
        encoder.load_state_dict(torch.load(encoder_file))
        decoder_file = os.path.join(args.load_folder, 'decoder.pt')
        decoder.load_state_dict(torch.load(decoder_file))

        args.save_folder = False

    # ===================================
    # set up training parameters
    # ===================================
    if args.optimizer == 'Adam':
        optimizer = optim.Adam(list(encoder.parameters()) + list(decoder.parameters()), lr=args.lr)
    elif args.optimizer == 'LBFGS':
        optimizer = optim.LBFGS(list(encoder.parameters()) + list(decoder.parameters()),
                                lr=args.lr)
    elif args.optimizer == 'SGD':
        optimizer = optim.SGD(list(encoder.parameters()) + list(decoder.parameters()),
                              lr=args.lr)

    scheduler = lr_scheduler.StepLR(optimizer, step_size=args.lr_decay,
                                    gamma=args.gamma)

    # Linear indices of an upper triangular mx, used for acc calculation
    triu_indices = get_triu_offdiag_indices(args.data_variable_size)
    tril_indices = get_tril_offdiag_indices(args.data_variable_size)

    if args.prior:
        prior = np.array([0.91, 0.03, 0.03, 0.03])  # hard coded for now
        print("Using prior")
        print(prior)
        log_prior = torch.DoubleTensor(np.log(prior))
        log_prior = torch.unsqueeze(log_prior, 0)
        log_prior = torch.unsqueeze(log_prior, 0)
        log_prior = Variable(log_prior)

        if args.cuda:
            log_prior = log_prior.cuda()

    if args.cuda:
        encoder.cuda()
        decoder.cuda()
        rel_rec = rel_rec.cuda()
        rel_send = rel_send.cuda()
        triu_indices = triu_indices.cuda()
        tril_indices = tril_indices.cuda()

    rel_rec = Variable(rel_rec)
    rel_send = Variable(rel_send)

    # ===================================
    # training:
    # ===================================
    def train(train_loader, epoch, best_val_loss, best_graph, ground_truth_G, lambda_A, c_A, optimizer):
        t = time.time()
        nll_train = []
        kl_train = []
        mse_train = []
        shd_trian = []

        encoder.train()
        decoder.train()
        scheduler.step()

        # update optimizer
        optimizer, lr = update_optimizer(optimizer, args.lr, c_A)

        for batch_idx, (data, relations) in enumerate(train_loader):

            if args.cuda:
                data, relations = data.cuda(), relations.cuda()
            data, relations = Variable(data).double(), Variable(relations).double()

            # reshape data
            relations = relations.unsqueeze(2)

            optimizer.zero_grad()

            enc_x, logits, origin_A, adj_A_tilt_encoder, z_gap, z_positive, myA, Wa = encoder(data, rel_rec,
                                                                                              rel_send)  # logits is of size: [num_sims, z_dims]
            edges = logits

            dec_x, output, adj_A_tilt_decoder = decoder(data, edges, args.data_variable_size * args.x_dims, rel_rec,
                                                        rel_send, origin_A, adj_A_tilt_encoder, Wa)

            if torch.sum(output != output):
                print('nan error\n')

            target = data
            preds = output
            variance = 0.

            # reconstruction accuracy loss
            loss_nll = nll_gaussian(preds, target, variance)

            # KL loss
            loss_kl = kl_gaussian_sem(logits)

            # ELBO loss:
            loss = loss_kl + loss_nll

            # add A loss
            one_adj_A = origin_A  # torch.mean(adj_A_tilt_decoder, dim =0)
            sparse_loss = args.tau_A * torch.sum(torch.abs(one_adj_A))

            # other loss term
            if args.use_A_connect_loss:
                connect_gap = A_connect_loss(one_adj_A, args.graph_threshold, z_gap)
                loss += lambda_A * connect_gap + 0.5 * c_A * connect_gap * connect_gap

            if args.use_A_positiver_loss:
                positive_gap = A_positive_loss(one_adj_A, z_positive)
                loss += .1 * (lambda_A * positive_gap + 0.5 * c_A * positive_gap * positive_gap)

            # compute h(A)
            h_A = _h_A(origin_A, args.data_variable_size)
            loss += lambda_A * h_A + 0.5 * c_A * h_A * h_A + 100. * torch.trace(
                origin_A * origin_A) + sparse_loss  # +  0.01 * torch.sum(variance * variance)

            loss.backward()
            loss = optimizer.step()

            myA.data = stau(myA.data, args.tau_A * lr)

            if torch.sum(origin_A != origin_A):
                print('nan error\n')

            # compute metrics
            graph = origin_A.data.clone().numpy()
            graph[np.abs(graph) < args.graph_threshold] = 0

            fdr, tpr, fpr, shd, nnz = count_accuracy(ground_truth_G, nx.DiGraph(graph))

            mse_train.append(F.mse_loss(preds, target).item())
            nll_train.append(loss_nll.item())
            kl_train.append(loss_kl.item())
            shd_trian.append(shd)

        if (np.mean(kl_train) + np.mean(nll_train)) < best_val_loss:
            best_graph = np.copy(graph)

        print(h_A.item())
        nll_val = []
        acc_val = []
        kl_val = []
        mse_val = []

        print('Epoch: {:04d}'.format(epoch),
              'nll_train: {:.10f}'.format(np.mean(nll_train)),
              'kl_train: {:.10f}'.format(np.mean(kl_train)),
              'ELBO_loss: {:.10f}'.format(np.mean(kl_train) + np.mean(nll_train)),
              'mse_train: {:.10f}'.format(np.mean(mse_train)),
              'shd_trian: {:.10f}'.format(np.mean(shd_trian)),
              'time: {:.4f}s'.format(time.time() - t))
        if args.save_folder and np.mean(nll_val) < best_val_loss:
            torch.save(encoder.state_dict(), encoder_file)
            torch.save(decoder.state_dict(), decoder_file)
            print('Best model so far, saving...')
            print('Epoch: {:04d}'.format(epoch),
                  'nll_train: {:.10f}'.format(np.mean(nll_train)),
                  'kl_train: {:.10f}'.format(np.mean(kl_train)),
                  'ELBO_loss: {:.10f}'.format(np.mean(kl_train) + np.mean(nll_train)),
                  'mse_train: {:.10f}'.format(np.mean(mse_train)),
                  'shd_trian: {:.10f}'.format(np.mean(shd_trian)),
                  'time: {:.4f}s'.format(time.time() - t), file=log)
            log.flush()

        if 'graph' not in vars():
            print('error on assign')

        return np.mean(np.mean(kl_train) + np.mean(nll_train)), np.mean(nll_train), np.mean(mse_train), graph, origin_A, best_graph

    # ===================================
    # main
    # ===================================

    t_total = time.time()
    best_ELBO_loss = np.inf
    best_NLL_loss = np.inf
    best_MSE_loss = np.inf
    best_epoch = 0
    best_ELBO_graph = []
    best_NLL_graph = []
    best_MSE_graph = []
    best_graph = np.ones((num_nodes, num_nodes)) - np.eye(num_nodes)
    # optimizer step on hyparameters
    c_A = args.c_A
    lambda_A = args.lambda_A
    h_A_new = torch.tensor(1.)
    h_tol = args.h_tol
    k_max_iter = int(args.k_max_iter)
    h_A_old = np.inf

    try:
        flag_max_iter = True
        for step_k in range(k_max_iter):
            while c_A < 1e+20:
                for epoch in range(args.epochs):
                    ELBO_loss, NLL_loss, MSE_loss, graph, origin_A, best_graph = train(train_loader, epoch,
                                                                                       best_ELBO_loss, best_graph,
                                                                                       ground_truth_G, lambda_A, c_A,
                                                                                       optimizer)
                    if ELBO_loss < best_ELBO_loss:
                        best_ELBO_loss = ELBO_loss
                        best_epoch = epoch
                        best_ELBO_graph = graph

                    if NLL_loss < best_NLL_loss:
                        best_NLL_loss = NLL_loss
                        best_epoch = epoch
                        best_NLL_graph = graph

                    if MSE_loss < best_MSE_loss:
                        best_MSE_loss = MSE_loss
                        best_epoch = epoch
                        best_MSE_graph = graph

                print("Optimization Finished!")
                print("Best Epoch: {:04d}".format(best_epoch))
                if ELBO_loss > 2 * best_ELBO_loss:
                    break

                # update parameters
                A_new = origin_A.data.clone()
                h_A_new = _h_A(A_new, args.data_variable_size)
                if h_A_new.item() > 0.25 * h_A_old:
                    c_A *= 10
                else:
                    break

                # update parameters
                # h_A, adj_A are computed in loss anyway, so no need to store
            h_A_old = h_A_new.item()
            lambda_A += c_A * h_A_new.item()

            if h_A_new.item() <= h_tol:
                flag_max_iter = False
                break

        if args.save_folder:
            print("Best Epoch: {:04d}".format(best_epoch), file=log)
            log.flush()

        return best_graph, rel_rec, rel_send, encoder, decoder, train_loader, flag_max_iter

    except KeyboardInterrupt:
        # print the best anway
        print(best_ELBO_graph)
        print(nx.to_numpy_array(ground_truth_G))
        fdr, tpr, fpr, shd, nnz = count_accuracy(ground_truth_G, nx.DiGraph(best_ELBO_graph))
        print('Best ELBO Graph Accuracy: fdr', fdr, ' tpr ', tpr, ' fpr ', fpr, 'shd', shd, 'nnz', nnz)

        print(best_NLL_graph)
        print(nx.to_numpy_array(ground_truth_G))
        fdr, tpr, fpr, shd, nnz = count_accuracy(ground_truth_G, nx.DiGraph(best_NLL_graph))
        print('Best NLL Graph Accuracy: fdr', fdr, ' tpr ', tpr, ' fpr ', fpr, 'shd', shd, 'nnz', nnz)

        print(best_MSE_graph)
        print(nx.to_numpy_array(ground_truth_G))
        fdr, tpr, fpr, shd, nnz = count_accuracy(ground_truth_G, nx.DiGraph(best_MSE_graph))
        print('Best MSE Graph Accuracy: fdr', fdr, ' tpr ', tpr, ' fpr ', fpr, 'shd', shd, 'nnz', nnz)


def retrain(args, train_data, valid_data, gt_adjacency, mask_A):
    args.cuda = not args.no_cuda and torch.cuda.is_available()
    args.factor = not args.no_factor
    args.save_folder = args.exp_path

    if args.dynamic_graph:
        print("Testing with dynamically re-computed graph.")

    # Save model and meta-data. Always saves in a new sub-folder.
    if args.save_folder:
        exp_counter = 0
        now = datetime.datetime.now()
        timestamp = now.isoformat()
        save_folder = '{}/exp{}/'.format(args.save_folder, timestamp)
        # safe_name = save_folder.text.replace('/', '_')
        os.makedirs(save_folder)
        meta_file = os.path.join(save_folder, 'metadata.pkl')
        encoder_file = os.path.join(save_folder, 'encoder.pt')
        decoder_file = os.path.join(save_folder, 'decoder.pt')

        log_file = os.path.join(save_folder, 'log.txt')
        log = open(log_file, 'w')

        pickle.dump({'args': args}, open(meta_file, "wb"))
    else:
        print("WARNING: No save_folder provided!" +
              "Testing (within this script) will throw an error.")

    # ================================================
    # get data: experiments = {synthetic SEM, ALARM}
    # ================================================
    train_loader, ground_truth_G = load_numpy_data(args, train_data, gt_adjacency)
    if valid_data is not None:
        valid_loader, _ = load_numpy_data(args, valid_data, gt_adjacency)
    else:
        valid_loader = None

    # ===================================
    # load modules
    # ===================================
    # Generate off-diagonal interaction graph
    off_diag = np.ones([args.data_variable_size, args.data_variable_size]) - np.eye(args.data_variable_size)

    rel_rec = np.array(encode_onehot(np.where(off_diag)[1]), dtype=np.float64)
    rel_send = np.array(encode_onehot(np.where(off_diag)[0]), dtype=np.float64)
    rel_rec = torch.DoubleTensor(rel_rec)
    rel_send = torch.DoubleTensor(rel_send)

    # add adjacency matrix A
    num_nodes = args.data_variable_size
    adj_A = np.zeros((num_nodes, num_nodes))

    if args.encoder == 'mlp':
        encoder = MLPEncoder(args.data_variable_size * args.x_dims, args.x_dims, args.encoder_hidden,
                             int(args.z_dims), adj_A, mask_A,
                             batch_size=args.batch_size,
                             do_prob=args.encoder_dropout, factor=args.factor).double()
    elif args.encoder == 'sem':
        encoder = SEMEncoder(args.data_variable_size * args.x_dims, args.encoder_hidden,
                             int(args.z_dims), adj_A,
                             batch_size=args.batch_size,
                             do_prob=args.encoder_dropout, factor=args.factor).double()

    if args.decoder == 'mlp':
        decoder = MLPDecoder(args.data_variable_size * args.x_dims,
                             args.z_dims, args.x_dims, encoder,
                             data_variable_size=args.data_variable_size,
                             batch_size=args.batch_size,
                             n_hid=args.decoder_hidden,
                             do_prob=args.decoder_dropout).double()
    elif args.decoder == 'sem':
        decoder = SEMDecoder(args.data_variable_size * args.x_dims,
                             args.z_dims, 2, encoder,
                             data_variable_size=args.data_variable_size,
                             batch_size=args.batch_size,
                             n_hid=args.decoder_hidden,
                             do_prob=args.decoder_dropout).double()

    if args.load_folder:
        encoder_file = os.path.join(args.load_folder, 'encoder.pt')
        encoder.load_state_dict(torch.load(encoder_file))
        decoder_file = os.path.join(args.load_folder, 'decoder.pt')
        decoder.load_state_dict(torch.load(decoder_file))

        args.save_folder = False

    # ===================================
    # set up training parameters
    # ===================================
    if args.optimizer == 'Adam':
        optimizer = optim.Adam(list(encoder.parameters()) + list(decoder.parameters()), lr=args.lr)
    elif args.optimizer == 'LBFGS':
        optimizer = optim.LBFGS(list(encoder.parameters()) + list(decoder.parameters()),
                                lr=args.lr)
    elif args.optimizer == 'SGD':
        optimizer = optim.SGD(list(encoder.parameters()) + list(decoder.parameters()),
                              lr=args.lr)

    scheduler = lr_scheduler.StepLR(optimizer, step_size=args.lr_decay,
                                    gamma=args.gamma)

    # Linear indices of an upper triangular mx, used for acc calculation
    triu_indices = get_triu_offdiag_indices(args.data_variable_size)
    tril_indices = get_tril_offdiag_indices(args.data_variable_size)

    if args.prior:
        prior = np.array([0.91, 0.03, 0.03, 0.03])  # hard coded for now
        print("Using prior")
        print(prior)
        log_prior = torch.DoubleTensor(np.log(prior))
        log_prior = torch.unsqueeze(log_prior, 0)
        log_prior = torch.unsqueeze(log_prior, 0)
        log_prior = Variable(log_prior)

        if args.cuda:
            log_prior = log_prior.cuda()

    if args.cuda:
        encoder.cuda()
        decoder.cuda()
        rel_rec = rel_rec.cuda()
        rel_send = rel_send.cuda()
        triu_indices = triu_indices.cuda()
        tril_indices = tril_indices.cuda()

    rel_rec = Variable(rel_rec)
    rel_send = Variable(rel_send)

    # ===================================
    # training:
    # ===================================
    def train(train_loader, optimizer):
        t = time.time()
        nll_train = []
        kl_train = []
        mse_train = []
        shd_trian = []

        encoder.train()
        decoder.train()
        scheduler.step()

        # update optimizer
        optimizer, lr = update_optimizer(optimizer, args.lr, 1.)

        for batch_idx, (data, relations) in enumerate(train_loader):

            if args.cuda:
                data, relations = data.cuda(), relations.cuda()
            data, relations = Variable(data).double(), Variable(relations).double()

            # reshape data
            relations = relations.unsqueeze(2)

            optimizer.zero_grad()

            enc_x, logits, origin_A, adj_A_tilt_encoder, z_gap, z_positive, myA, Wa = encoder(data, rel_rec,
                                                                                              rel_send)  # logits is of size: [num_sims, z_dims]
            edges = logits

            dec_x, output, adj_A_tilt_decoder = decoder(data, edges, args.data_variable_size * args.x_dims, rel_rec,
                                                        rel_send, origin_A, adj_A_tilt_encoder, Wa)

            if torch.sum(output != output):
                print('nan error\n')

            target = data
            preds = output
            variance = 0.

            # reconstruction accuracy loss
            loss_nll = nll_gaussian(preds, target, variance)

            # KL loss
            loss_kl = kl_gaussian_sem(logits)

            # ELBO loss:
            loss = loss_kl + loss_nll

            # add A loss
            one_adj_A = origin_A  # torch.mean(adj_A_tilt_decoder, dim =0)
            sparse_loss = args.tau_A * torch.sum(torch.abs(one_adj_A))

            loss.backward()
            optimizer.step()

            myA.data = stau(myA.data, args.tau_A * lr)

            if torch.sum(origin_A != origin_A):
                print('nan error\n')

        return loss

    # ===================================
    # main
    # ===================================

    t_total = time.time()
    # optimizer step on hyparameters
    k_max_iter = int(args.k_max_iter)
    best_valid_score = -np.inf
    best_train_score = None
    full_patience = 5
    patience = full_patience

    flag_max_iter = True
    for epoch in range(100*args.epochs):
        train_score = -train(train_loader, optimizer)
        if valid_loader is not None:
            valid_score = evaluate_score(encoder, decoder, valid_loader, args.graph_threshold, rel_rec, rel_send,
                                         num_nodes, args)
        else:
            valid_score = train_score

        if valid_score > best_valid_score:
            best_valid_score = valid_score
            best_train_score = evaluate_score(encoder, decoder, train_loader, args.graph_threshold, rel_rec, rel_send,
                                              num_nodes, args)
            patience = full_patience
        else:
            patience -= 1
        if patience == 0:
            flag_max_iter = False
            break

    if valid_loader is None:
        best_valid_score = None

    return best_train_score, best_valid_score, flag_max_iter

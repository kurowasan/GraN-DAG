"""The full NOTEARS algorithm with l1 regularization.

Defines two functions:
- notears(): the vanilla python implementation of the algorithm
- notears_live(): same algorithm with live plotting in jupyter notebook

The full NOTEARS algorithm works with both large n and small n.
"""
import numpy as np
import networkx as nx
from bokeh.plotting import figure
from bokeh.models import ColumnDataSource, LinearColorMapper
from bokeh.transform import transform
from bokeh.palettes import RdBu11 as Palette
from bokeh.models import HoverTool
from bokeh.layouts import gridplot
from bokeh.io import show, push_notebook

import cppext
import torch

def notears(X: np.ndarray,
            pns_mask: np.ndarray,
            lambda1: float,
            max_iter: int = 100,
            h_tol: float = 1e-8,
            w_threshold: float = 0.3) -> np.ndarray:
    """Solve min_W F(W; X) s.t. h(W) = 0 using augmented Lagrangian.

    Args:
        X: [n,d] sample matrix
        lambda1: l1 regularization parameter
        max_iter: max number of dual ascent steps
        h_tol: exit if |h(w)| <= h_tol
        w_threshold: fixed threshold for edge weights

    Returns:
        W_est: [d,d] estimate
    """
    n, d = X.shape
    w_est, w_new = np.zeros(d * d), np.zeros(d * d)
    pns_mask = pns_mask.astype(np.float_)
    rho, alpha, h, h_new = 1.0, 0.0, np.inf, np.inf
    flag_max_iter = True
    for i in range(max_iter):
        while rho < 1e+20:
            w_new = cppext.minimize_subproblem(w_est, X, pns_mask, rho, alpha, lambda1)
            h_new = cppext.h_func(w_new)
            if h_new > 0.25 * h:
                rho *= 10
            else:
                break
        w_est, h = w_new, h_new
        alpha += rho * h
        if h <= h_tol:
            flag_max_iter = False
            break
    w_est[np.abs(w_est) < w_threshold] = 0
    return w_est.reshape([d, d]), flag_max_iter


def notears_live(G: nx.DiGraph,
                 X: np.ndarray,
                 pns_mask: np.ndarray,
                 lambda1: float,
                 max_iter: int = 100,
                 h_tol: float = 1e-8,
                 w_threshold: float = 0.3) -> np.ndarray:
    """Monitor the optimization progress live in notebook.

    Args:
        G: ground truth graph
        X: [n,d] sample matrix
        lambda1: l1 regularization parameter
        max_iter: max number of dual ascent steps
        h_tol: exit if |h(w)| <= h_tol
        w_threshold: fixed threshold for edge weights

    Returns:
        W_est: [d,d] estimate
    """
    # initialization
    n, d = X.shape
    w_est, w_new = np.zeros(d * d), np.zeros(d * d)
    rho, alpha, h, h_new = 1.0, 0.0, np.inf, np.inf

    # ground truth
    w_true = nx.to_numpy_array(G).flatten()

    # progress, stream
    progress_data = {key:[] for key in ['step', 'F', 'h',
                                        'rho', 'alpha', 'l2_dist']}
    progress_source = ColumnDataSource(data=progress_data)

    # heatmap, patch
    ids = [str(i) for i in range(d)]
    all_ids = np.tile(ids, [d, 1])
    row = all_ids.T.flatten()
    col = all_ids.flatten()
    heatmap_data = {'row': row, 'col': col,
                    'w_true': w_true, 'w_est': w_est, 'w_diff': w_true - w_est}
    heatmap_source = ColumnDataSource(data=heatmap_data)
    mapper = LinearColorMapper(palette=Palette, low=-2, high=2)

    # common tools
    tools = 'crosshair,save,reset'

    # F(w_est) vs step
    F_true = cppext.F_func(w_true, X, lambda1)
    fig0 = figure(plot_width=270, plot_height=240,
                  y_axis_type='log', tools=tools)
    fig0.ray(0, F_true, length=0, angle=0, color='green',
             line_dash='dashed', line_width=2, legend='F(w_true)')
    fig0.line('step', 'F', source=progress_source,
              color='red', line_width=2, legend='F(w_est)')
    fig0.title.text = "Objective"
    fig0.xaxis.axis_label = "step"
    fig0.legend.location = "bottom_left"
    fig0.legend.background_fill_alpha = 0.5
    fig0.add_tools(HoverTool(tooltips=[("step", "@step"),
                                       ("F", "@F"),
                                       ("F_true", '%.6g' % F_true)],
                             mode='vline'))

    # h(w_est) vs step
    fig1 = figure(plot_width=280, plot_height=240,
                  y_axis_type='log', tools=tools)
    fig1.line('step', 'h', source=progress_source,
              color='magenta', line_width=2, legend='h(w_est)')
    fig1.title.text = "Constraint"
    fig1.xaxis.axis_label = "step"
    fig1.legend.location = "bottom_left"
    fig1.legend.background_fill_alpha = 0.5
    fig1.add_tools(HoverTool(tooltips=[("step", "@step"),
                                       ("h", "@h"),
                                       ("rho", "@rho"),
                                       ("alpha", "@alpha")],
                             mode='vline'))

    # ||w_true - w_est|| vs step
    fig2 = figure(plot_width=270, plot_height=240,
                  y_axis_type='log', tools=tools)
    fig2.line('step', 'l2_dist', source=progress_source,
              color='blue', line_width=2)
    fig2.title.text = "L2 distance to W_true"
    fig2.xaxis.axis_label = "step"
    fig2.add_tools(HoverTool(tooltips=[("step", "@step"),
                                       ("w_est", "@l2_dist")],
                             mode='vline'))

    # heatmap of w_true
    fig3 = figure(plot_width=270, plot_height=240,
                  x_range=ids, y_range=list(reversed(ids)), tools=tools)
    fig3.rect(x='col', y='row', width=1, height=1, source=heatmap_source,
              line_color=None, fill_color=transform('w_true', mapper))
    fig3.title.text = 'W_true'
    fig3.axis.visible = False
    fig3.add_tools(HoverTool(tooltips=[("row, col", "@row, @col"),
                                       ("w_true", "@w_true")]))

    # heatmap of w_est
    fig4 = figure(plot_width=280, plot_height=240,
                  x_range=ids, y_range=list(reversed(ids)), tools=tools)
    fig4.rect(x='col', y='row', width=1, height=1, source=heatmap_source,
              line_color=None, fill_color=transform('w_est', mapper))
    fig4.title.text = 'W_est'
    fig4.axis.visible = False
    fig4.add_tools(HoverTool(tooltips=[("row, col", "@row, @col"),
                                       ("w_est", "@w_est")]))

    # heatmap of w_true - w_est
    fig5 = figure(plot_width=270, plot_height=240,
                  x_range=ids, y_range=list(reversed(ids)), tools=tools)
    fig5.rect(x='col', y='row', width=1, height=1, source=heatmap_source,
               line_color=None, fill_color=transform('w_diff', mapper))
    fig5.title.text = 'W_true - W_est'
    fig5.axis.visible = False
    fig5.add_tools(HoverTool(tooltips=[("row, col", "@row, @col"),
                                       ("w_diff", "@w_diff")]))

    # display figures as grid
    grid = gridplot([[fig0, fig1, fig2],
                     [fig3, fig4, fig5]], merge_tools=False)
    handle = show(grid, notebook_handle=True)

    # enter main loop
    for it in range(max_iter):
        while rho < 1e+20:
            w_new = cppext.minimize_subproblem(w_est, X, pns_mask, rho, alpha, lambda1)
            h_new = cppext.h_func(w_new)
            if h_new > 0.25 * h:
                rho *= 10
            else:
                break
        w_est, h = w_new, h_new
        alpha += rho * h
        # update figures
        progress_source.stream({'step': [it],
                                'F': [cppext.F_func(w_est, X, lambda1)],
                                'h': [h],
                                'rho': [rho],
                                'alpha': [alpha],
                                'l2_dist': [np.linalg.norm(w_est - w_true)],})
        heatmap_source.patch({'w_est': [(slice(d * d), w_est)],
                              'w_diff': [(slice(d * d), w_true - w_est)]})
        push_notebook(handle=handle)
        # check termination of main loop
        if h <= h_tol:
            break

    # final threshold
    w_est[np.abs(w_est) < w_threshold] = 0
    return w_est.reshape([d, d])

class LinearMasked(torch.nn.Module):
    def __init__(self, adj):
        super(LinearMasked, self).__init__()
        self.weights = torch.nn.Parameter(torch.Tensor(np.zeros_like(adj)))
        self.adj = torch.Tensor(adj)

    def forward(self, data):
        out = torch.matmul(data, self.adj * self.weights)
        return 0.5 * torch.sum((data - out)**2) / data.size(0)

def retrain(adj, train_data, valid_data, max_iter):

    model = LinearMasked(adj)
    optimizer = torch.optim.Adagrad(model.parameters())
    best_valid_score = -np.inf
    full_patience = 100
    flag_max_iter = True

    for iter in range(max_iter):
        x, _ = train_data.sample(train_data.num_samples)

        loss = model(x)
        train_score = -loss.item()

        optimizer.zero_grad()
        loss.backward()
        optimizer.step()

        if valid_data is not None:
            with torch.no_grad():
                x, _ = valid_data.sample(valid_data.num_samples)
                valid_score = -model(x).item()
        else:
            valid_score = train_score

        if valid_score > best_valid_score + 1e-4:
            best_valid_score = valid_score
            best_train_score = train_score
            patience = full_patience
        else:
            patience -= 1

        if iter % 100 == 0:
            print("Iteration: {}, score_train: {} , score_valid : {}, best_valid_score {}, patience: {}".format(iter, train_score, valid_score, best_valid_score, patience))
        if patience == 0:
            flag_max_iter = False
            break

    return best_train_score, best_valid_score, flag_max_iter



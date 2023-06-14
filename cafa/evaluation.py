import numpy as np
import logging
import pandas as pd
import multiprocessing as mp


# Computes the root terms in the dag
def get_roots_idx(dag):
    return np.where(dag.sum(axis=1) == 0)[0]


# Computes the leaf terms in the dag
def get_leafs_idx(dag):
    return np.where(dag.sum(axis=0) == 0)[0]


# Return a mask for all the predictions (matrix) >= tau
def solidify_prediction(pred, tau):
    return pred >= tau


# computes the f metric for each precision and recall in the input arrays
def compute_f(pr, rc):
    n = 2 * pr * rc
    d = pr + rc
    return np.divide(n, d, out=np.zeros_like(n, dtype=float), where=d != 0)


def compute_s(ru, mi):
    return np.sqrt(ru**2 + mi**2)
    # return np.where(np.isnan(ru), mi, np.sqrt(ru + np.nan_to_num(mi)))


def compute_metrics_(tau_arr, g, pred, toi, n_gt, wn_gt=None, ic_arr=None):

    metrics = np.zeros((len(tau_arr), 7), dtype='float')  # cov, pr, rc, wpr, wrc, ru, mi

    for i, tau in enumerate(tau_arr):

        p = solidify_prediction(pred.matrix[:, toi], tau)

        # number of proteins with at least one term predicted with score >= tau
        metrics[i, 0] = (p.sum(axis=1) > 0).sum()

        # Terms subsets
        intersection = np.logical_and(p, g)  # TP

        # Subsets size
        n_pred = p.sum(axis=1)
        n_intersection = intersection.sum(axis=1)

        # Precision, recall
        metrics[i, 1] = np.divide(n_intersection, n_pred, out=np.zeros_like(n_intersection, dtype='float'),
                                  where=n_pred > 0).sum()
        metrics[i, 2] = np.divide(n_intersection, n_gt, out=np.zeros_like(n_gt, dtype='float'), where=n_gt > 0).sum()

        if ic_arr is not None:
            # Terms subsets
            remaining = np.logical_and(np.logical_not(p), g)  # FN --> not predicted but in the ground truth
            mis = np.logical_and(p, np.logical_not(g))  # FP --> predicted but not in the ground truth

            # Weighted precision, recall
            wn_pred = (p * ic_arr[toi]).sum(axis=1)
            wn_intersection = (intersection * ic_arr[toi]).sum(axis=1)

            metrics[i, 3] = np.divide(wn_intersection, wn_pred, out=np.zeros_like(n_intersection, dtype='float'),
                                      where=n_pred > 0).sum()
            metrics[i, 4] = np.divide(wn_intersection, wn_gt, out=np.zeros_like(n_intersection, dtype='float'),
                                      where=n_gt > 0).sum()

            # Misinformation, remaining uncertainty
            metrics[i, 5] = (remaining * ic_arr[toi]).sum(axis=1).sum()
            metrics[i, 6] = (mis * ic_arr[toi]).sum(axis=1).sum()
    return metrics


def compute_metrics(pred, gt, toi, tau_arr, ic_arr=None, n_cpu=0):
    """
    Takes the prediction and the ground truth and for each threshold in tau_arr
    calculates the confusion matrix and returns the coverage,
    precision, recall, remaining uncertainty and misinformation.
    Toi is the list of terms (indexes) to be considered
    """
    g = gt.matrix[:, toi]
    n_gt = g.sum(axis=1)
    wn_gt = None
    if ic_arr is not None:
        wn_gt = (g * ic_arr[toi]).sum(axis=1)

    # Parallelization
    if n_cpu == 0:
        n_cpu = mp.cpu_count()

    arg_lists = [[tau_arr, g, pred, toi, n_gt, wn_gt, ic_arr] for tau_arr in np.array_split(tau_arr, n_cpu)]
    with mp.Pool(processes=n_cpu) as pool:
        metrics = np.concatenate(pool.starmap(compute_metrics_, arg_lists), axis=0)

    return pd.DataFrame(metrics, columns=["cov", "pr", "rc", "wpr", "wrc", "ru", "mi"])


def evaluate_prediction(prediction, gt, ontologies, tau_arr, normalization='cafa', n_cpu=0):

    dfs = []
    for p in prediction:
        ns = p.namespace
        ne = np.full(len(tau_arr), gt[ns].matrix.shape[0])

        ont = [o for o in ontologies if o.namespace == ns][0]

        # cov, pr, rc, wpr, wrc, ru, mi
        metrics = compute_metrics(p, gt[ns], ont.toi, tau_arr, ont.ia, n_cpu)

        for column in ["pr", "rc", "wpr", "wrc", "ru", "mi"]:
            if normalization == 'gt' or (column in ["rc", "wrc"] and normalization == 'cafa'):
                metrics[column] = np.divide(metrics[column], ne, out=np.zeros_like(metrics[column], dtype='float'), where=ne > 0)
            else:
                metrics[column] = np.divide(metrics[column], metrics["cov"], out=np.zeros_like(metrics[column], dtype='float'), where=metrics["cov"] > 0)

        metrics['ns'] = [ns] * len(tau_arr)
        metrics['tau'] = tau_arr
        metrics['cov'] = np.divide(metrics['cov'], ne, out=np.zeros_like(metrics['cov'], dtype='float'), where=ne > 0)
        metrics['f'] = compute_f(metrics['pr'], metrics['rc'])
        metrics['wf'] = compute_f(metrics['wpr'], metrics['wrc'])
        metrics['s'] = compute_s(metrics['ru'], metrics['mi'])

        dfs.append(metrics)

    return pd.concat(dfs)

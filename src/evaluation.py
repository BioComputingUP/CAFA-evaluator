import numpy as np
import logging
import pandas as pd


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


def compute_metrics(pred, gt, toi, tau_arr, ic_arr=None):
    """
    takes the prediction, the ground truth and returns the sum of the terms in order to compute precision and recall
    In the case of precision a tuple (precison_sum, m_tau),
    m_tau is the number of rows with at least one score >= tau
    """
    pr_list = np.zeros(len(tau_arr), dtype='float')
    rc_list = np.zeros(len(tau_arr), dtype='float')
    ru_list = np.zeros(len(tau_arr), dtype='float')
    mi_list = np.zeros(len(tau_arr), dtype='float')
    wpr_list = np.zeros(len(tau_arr), dtype='float')
    wrc_list = np.zeros(len(tau_arr), dtype='float')
    cov_list = np.zeros(len(tau_arr), dtype='int')

    g = gt.matrix[:, toi]
    n_gt = g.sum(axis=1)
    if ic_arr is not None:
        wn_gt = (g * ic_arr[toi]).sum(axis=1)

    for i, tau in enumerate(tau_arr):
        p = solidify_prediction(pred.matrix[:, toi], tau)
        cov_list[i] = (p.sum(axis=1) > 0).sum()  # number of proteins with at least one term predicted with score >= tau

        # Terms subsets
        intersection = np.logical_and(p, g)  # TP
        remaining = np.logical_and(np.logical_not(p), g)  # FN --> not predicted but in the ground truth
        mis = np.logical_and(p, np.logical_not(g))  # FP --> predicted but not in the ground truth

        # Subsets size
        n_pred = p.sum(axis=1)
        n_intersection = intersection.sum(axis=1)

        # Precision, recall
        pr_list[i] = np.divide(n_intersection, n_pred, out=np.zeros_like(n_intersection, dtype='float'), where=n_pred > 0).sum()
        rc_list[i] = np.divide(n_intersection, n_gt, out=np.zeros_like(n_gt, dtype='float'), where=n_gt > 0).sum()

        if ic_arr is not None:

            # Weighted precision, recall
            wn_pred = (p * ic_arr[toi]).sum(axis=1)
            wn_intersection = (intersection * ic_arr[toi]).sum(axis=1)

            wpr_list[i] = np.divide(wn_intersection, wn_pred, out=np.zeros_like(n_intersection, dtype='float'), where=n_pred > 0).sum()
            wrc_list[i] = np.divide(wn_intersection, wn_gt, out=np.zeros_like(n_intersection, dtype='float'), where=n_gt > 0).sum()

            # Misinformation, remaining uncertainty
            ru_list[i] = (remaining * ic_arr[toi]).sum(axis=1).sum()
            mi_list[i] = (mis * ic_arr[toi]).sum(axis=1).sum()

    return [cov_list, pr_list, rc_list, wpr_list, wrc_list, ru_list, mi_list]


def evaluate_prediction(prediction, gt, ontologies, tau_arr, normalization='cafa'):

    dfs = []
    for p in prediction:
        ns = p.namespace
        ne = np.full(len(tau_arr), gt[ns].matrix.shape[0])

        ont = [o for o in ontologies if o.namespace == ns][0]

        # cov_list, pr_list, rc_list, wpr_list, wrc_list, ru_list, mi_list
        metrics = compute_metrics(p, gt[ns], ont.toi, tau_arr, ont.ia)
        cov_list = metrics[0]
        metrics = metrics[1:]

        # pr_list, rc_list, wpr_list, wrc_list, ru_list, mi_list
        # rc_list is treated differently in cafa normalization
        for i, metric in enumerate(metrics):
            if normalization == 'gt' or (i in [1, 3] and normalization == 'cafa'):
                metrics[i] = np.divide(metric, ne, out=np.zeros_like(metric, dtype='float'), where=ne > 0)
            else:  # pred and cafa
                metrics[i] = np.divide(metric, cov_list, out=np.zeros_like(metric, dtype='float'), where=cov_list > 0)
        pr_list, rc_list, wpr_list, wrc_list, ru_list, mi_list = metrics

        dfs.append(pd.DataFrame({'ns': [ns] * len(tau_arr),
                                 'tau': tau_arr,
                                 'cov': np.divide(cov_list, ne, out=np.zeros_like(cov_list, dtype='float'), where=ne > 0),
                                 'pr': pr_list,
                                 'rc': rc_list,
                                 'f': compute_f(pr_list, rc_list),
                                 'wpr': wpr_list,
                                 'wrc': wrc_list,
                                 'wf': compute_f(wpr_list, wrc_list),
                                 'mi': mi_list,
                                 'ru': ru_list,
                                 's':  compute_s(ru_list, mi_list)}))

    return pd.concat(dfs)

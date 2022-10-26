import numpy as np
import logging
import pandas as pd

from graph import propagate

# Tau array, used to compute metrics at different score thresholds
tau_arr = np.arange(0.01, 1, 0.01)


# Computes the root terms in the dag
def get_roots_idx(dag):
    return np.where(dag.sum(axis=1) == 0)[0]


# Computes the leaf terms in the dag
def get_leafs_idx(dag):
    return np.where(dag.sum(axis=0) == 0)[0]


# Return a mask for all the predictions (matrix) >= tau
def solidify_prediction(pred, tau):
    return pred >= tau


# Return the rows of the ground truth matrix that match the rows of the prediction matrix
def build_gt(preds, gt, gt_ids):
    g = np.zeros((len(preds.ids), gt.shape[1]), dtype='bool')
    for p_id, idx in preds.ids.items():
        if p_id in gt_ids:
            gt_idx = gt_ids[p_id]
            g[idx, :] = gt[gt_idx, :]
        else:
            logging.info("protein id:", p_id)
            raise Exception("A predicted protein id is not contained in the ground truth") 
    return g


# computes the f metric for each precision and recall in the input arrays
def compute_f(pr, rc):
    n = 2 * pr * rc
    d = pr + rc
    return np.divide(n, d, out=np.zeros_like(n, dtype=float), where=d != 0)


def compute_s(ru, mi):
    return np.sqrt(ru**2 + mi**2)
    # return np.where(np.isnan(ru), mi, np.sqrt(ru + np.nan_to_num(mi)))


def compute_metrics(pred, gt, tau_arr, toi, ic_arr=None):
    """
    takes the prediction, the ground truth and returns the sum of the terms in order to compute precision and recall
    In the case of precision a tuple (precison_sum, m_tau),
    m_tau is the number of rows with at least one score >= tau
    """
    pr_list = np.zeros((len(tau_arr), 2), dtype='float')
    rc_list = np.zeros(len(tau_arr), dtype='float')

    ru_list = np.zeros((len(tau_arr), 2))
    mi_list = np.zeros((len(tau_arr), 2))

    wpr_list = np.zeros((len(tau_arr), 2), dtype='float')
    wrc_list = np.zeros(len(tau_arr), dtype='float')

    for i, tau in enumerate(tau_arr):
        p = solidify_prediction(pred.matrix[:, toi], tau)
        g = build_gt(pred, gt.matrix, gt.ids)[:, toi]
        m_tau = ((p.sum(axis=1)) >= 1).sum()

        # Terms subsets
        intersection = np.logical_and(p, g)  # TP
        remaining = np.logical_and(np.logical_not(p), g)  # FN --> not predicted but in the ground truth
        mis = np.logical_and(p, np.logical_not(g))  # FP --> predicted but not in the ground truth

        # Subsets size
        n_gt = g.sum(axis=1)
        n_pred = p.sum(axis=1)
        n_intersection = intersection.sum(axis=1)

        # Precision, recall
        pr = np.divide(n_intersection, n_pred, out=np.zeros_like(n_intersection, dtype='float'), where=n_pred != 0).sum()
        rc = np.divide(n_intersection, n_gt, out=np.zeros_like(n_intersection, dtype='float'), where=n_gt != 0).sum()

        pr_list[i, :] = [pr, m_tau]
        rc_list[i] = rc

        if ic_arr is not None:

            # Weighted precision, recall
            wn_gt = (g * ic_arr).sum(axis=1)
            wn_pred = (p * ic_arr).sum(axis=1)
            wn_intersection = (intersection * ic_arr).sum(axis=1)

            wpr = np.divide(wn_intersection, wn_pred, out=np.zeros_like(wn_intersection, dtype='float'), where=wn_pred != 0).sum()
            wrc = np.divide(wn_intersection, wn_gt, out=np.zeros_like(wn_intersection, dtype='float'), where=wn_gt != 0).sum()

            wpr_list[i, :] = [wpr, m_tau]
            wrc_list[i] = wrc

            # Misinformation, remaining uncertainty
            ru = (remaining * ic_arr).sum(axis=1).sum()
            mi = (mis * ic_arr).sum(axis=1).sum()

            ru_list[i, :] = [ru, m_tau]
            mi_list[i, :] = [mi, m_tau]

    return pr_list, rc_list, wpr_list, wrc_list, ru_list, mi_list


def evaluate_prediction(prediction, gt, ontologies, ne):

    dfs = []
    for p in prediction:
        ns = p.namespace
        if ns in gt:
            ont = [o for o in ontologies if o.namespace == ns][0]

            _ = propagate(p, ont, ont.order)

            pr_sum, rc_sum, wpr_sum, wrc_sum, ru_sum, mi_sum = compute_metrics(p, gt[ns], tau_arr, ont.toi, ont.ia)

            # Precision recall
            n = pr_sum[:, 0]
            d = pr_sum[:, 1]
            _pr = np.divide(n, d, out=np.zeros_like(n, dtype='float'), where=d != 0)
            _rc = rc_sum / ne[ns]
            cov = d / ne[ns]

            # Mi, ru
            n = ru_sum[:, 0]
            d = ru_sum[:, 1]
            _ru = np.divide(n, d, out=np.zeros_like(n, dtype='float'), where=d != 0)

            n = mi_sum[:, 0]
            d = mi_sum[:, 1]
            _mi = np.divide(n, d, out=np.zeros_like(n, dtype='float'), where=d != 0)

            # Weighted precision recall
            n = wpr_sum[:, 0]
            d = wpr_sum[:, 1]
            _wpr = np.divide(n, d, out=np.zeros_like(n, dtype='float'), where=d != 0)
            _wrc = wrc_sum / ne[ns]

            dfs.append(pd.DataFrame({'ns': [ns] * len(tau_arr),
                                     'tau': tau_arr,
                                     'cov': cov,
                                     'pr': _pr,
                                     'rc': _rc,
                                     'f': compute_f(_pr, _rc),
                                     'wpr': _wpr,
                                     'wrc': _wrc,
                                     'wf': compute_f(_wpr, _wrc),
                                     'mi': _mi,
                                     'ru': _ru,
                                     's':  compute_s(_ru, _mi)}))

    return pd.concat(dfs)

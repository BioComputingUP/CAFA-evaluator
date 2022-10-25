import numpy as np
import matplotlib.pyplot as plt
import logging


# Compute the indexes of the term of interests in the dag
def get_toi_idx(dag):
    return np.where(dag.sum(axis=1) > 0)[0]


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
        intersection = np.logical_and(p, g)
        remaining = np.logical_and(np.logical_not(p), g)
        mis = np.logical_and(p, np.logical_not(g))

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


# TODO manage baselines colors
def plot_pr_rc(df, groups, out_file):
    """
    Plot precision recall curves
    """
    cmap = plt.get_cmap('tab10')
    fig, ax = plt.subplots(figsize=(10, 10))
    for method, df_m in df.groupby(level='method', sort=False):
        best = df_m.loc[df_m['f'].idxmax()]
        color = cmap.colors[groups.get_loc(best.name[1]) % len(cmap.colors)]
        # TODO baseline color is hardcoded
        if 'blast' in best.get('label', '').lower():
            color = (0, 0, 1)  # blue
        elif 'naive' in best.get('label', '').lower():
            color = (1, 0, 0)  # red
        ax.plot(df_m['rc'], df_m['pr'], '--' if best.get('is_baseline') else '-',
                label="{} (F={:.2f},C={:.2f})".format(best.get('label', method), best['f'], best['cov_f']), color=color)
        plt.plot(best['rc'], best['pr'], 'o', color=color)

    ax.legend()
    # ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', prop={'size': 15})
    # ax.legend(bbox_to_anchor=(0, 1), loc='lower left', prop={'size': 15}, ncol=int(math.sqrt(len(labels)))-1)
    ax.set_ylabel("Precision") #, fontsize=18)
    ax.set_xlabel("Recall") #, fontsize=18)
    plt.xlim([0, 1])
    plt.ylim([0, 1])
    plt.savefig(out_file, bbox_inches='tight')
    plt.close('all')


def plot_mi_ru(df, groups, out_file):
    """
    Plot misinformation remaining uncertainty curves
    """
    cmap = plt.get_cmap('tab10')
    fig, ax = plt.subplots(figsize=(10, 10))
    # plot the curves
    for method, df_m in df.groupby(level='method', sort=False):
        best = df_m.loc[df_m['s'].idxmin()]
        color = cmap.colors[groups.get_loc(best.name[1]) % len(cmap.colors)]
        # TODO baseline color is hardcoded
        if 'blast' in best.get('label', '').lower():
            color = (0, 0, 1)  # blue
        elif 'naive' in best.get('label', '').lower():
            color = (1, 0, 0)  # red
        ax.plot(df_m['ru'], df_m['mi'], '--' if best.get('is_baseline') else '-',
                label="{} (S={:.2f},C={:.2f})".format(best.get('label', method), best['s'], best['cov_s']), color=color)
        plt.plot(best['ru'], best['mi'], 'o', color=color)

    ax.legend()
    ax.set_ylabel("Misinformation")  # , fontsize=18)
    ax.set_xlabel("Remaining uncertainty")  # , fontsize=18)
    plt.xlim([0, 5])
    plt.ylim(bottom=0)
    plt.savefig(out_file, bbox_inches='tight')
    plt.close('all')

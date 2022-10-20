import numpy as np
import matplotlib.pyplot as plt
import logging
import pandas as pd

# compute the indexes of the term of interests in the dag
def get_toi_idx(dag):
    return np.where(dag.sum(axis=1) > 0)[0]


# computes the root terms in the dag
def get_roots_idx(dag):
    return np.where(dag.sum(axis=1) == 0)[0]


# computes the leaf terms in the dag
def get_leafs_idx(dag):
    return np.where(dag.sum(axis=0) == 0)[0]


# takes the prediction and a tau value and return a 0 1 matrix containing all the scores >= tau
def solidify_prediction(pred, tau):#, p, order, start_pos):
    return pred >= tau


# takes the rows of the gt matrix that match the rows of the prediction matrix
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


# computes the sum part of the precision for a prediction
def compute_precision_sum(pred, gt):
    n = np.logical_and(pred, gt).sum(axis=1)
    d = pred.sum(axis=1)
    s = np.divide(n, d, out=np.zeros_like(n, dtype='float'), where=d != 0).sum()
    return s


# nt is the number of target, if we use full eval mode nt is the number of benchmark proteins
# meanwhile if we use the partial evaluation mode nt is the number of proteins wich we gave in
# input to the predictor
def compute_recall_sum(pred, gt):
    n = np.logical_and(pred, gt).sum(axis=1)
    d = gt.sum(axis=1)
    s = np.divide(n, d, out=np.zeros_like(n, dtype='float'), where=d != 0).sum()
    return s


# takes the prediction, the ground truth and returns the sum of the terms in order to compute precision and recall
# in the case of precision a tuple precison_sum and m_tau is returned where m_tau is the number of rows with atleast 
# one score >= tau
def compute_metrics(pred, gt, tau_arr, toi, ic_arr):
    pr_list = np.zeros((len(tau_arr), 2), dtype='float')
    rc_list = np.zeros(len(tau_arr), dtype='float')

    ru_list = np.zeros((len(tau_arr), 2))
    mi_list = np.zeros((len(tau_arr), 2))

    i = 0
    for tau in tau_arr:
        p = solidify_prediction(pred.matrix[:, toi], tau)
        g = build_gt(pred, gt.matrix, gt.ids)[:, toi]
        m_tau = ((p.sum(axis=1)) >= 1).sum()

        pr = compute_precision_sum(p, g)
        rc = compute_recall_sum(p, g)
        pr_list[i, :] = [pr, m_tau]
        rc_list[i] = rc

        ru = compute_remaining_uncertainty(p, g, ic_arr[toi])
        mi = compute_misinformation(p, g, ic_arr[toi])
        ru_list[i, :] = [ru, m_tau]
        mi_list[i, :] = [mi, m_tau]

        i += 1
    return pr_list, rc_list, ru_list, mi_list


# computes the f metric for each precision and recall in the input arraysg
def compute_f(pr, rc):
    n = 2 * pr * rc
    d = pr + rc
    return np.divide(n, d, out=np.zeros_like(n, dtype=float), where=d != 0)


def compute_remaining_uncertainty(pred, gt, ic_arr):
    # n = np.logical_and(np.logical_not(pred), gt)
    return (np.logical_and(np.logical_not(pred), gt) * ic_arr).sum(axis=1).sum()


def compute_misinformation(pred, gt, ic_arr):
    return (np.logical_and(pred, np.logical_not(gt)) * ic_arr).sum(axis=1).sum()


def compute_s(ru, mi):
    return np.sqrt(ru**2 + mi**2)
    # return np.where(np.isnan(ru), mi, np.sqrt(ru + np.nan_to_num(mi)))


# plots a graph for each different namespace containing the curves for each (or the best based on f) team's prediction
def plot_pr_rc(data, out_folder):

    cmap = plt.get_cmap('tab10')
    groups = data.index.get_level_values('group').unique()

    for ns, df_ns in data.groupby(level='ns'):
        fig, ax = plt.subplots(figsize=(10, 10))

        # identify the best method for each group based on f max
        best_methods = []
        for group, df_g in df_ns.groupby(level='group'):
            best_methods.append(df_g['f'].idxmax())

        # sort methods based on f
        df_best_methods = df_ns.loc[best_methods].sort_values(by=['is_baseline', 'f'], ascending=[True, False])

        # sort threshold
        df_best_methods = df_best_methods.set_index('tau', append=True).groupby(level='method', group_keys=False, sort=False).apply(lambda x: x.sort_index(level='tau', ascending=False))

        # save to file
        df_best_methods.to_csv("{}/prrc_{}.tsv".format(out_folder, ns), sep='\t', float_format="%.3f")

        # plot the curves
        for method, df_m in df_best_methods.groupby(level='method', sort=False):
            best = df_m.loc[df_m['f'].idxmax()]
            color = cmap.colors[groups.get_loc(best.name[1]) % len(cmap.colors)]
            if 'blast' in best['label'].lower():
                color = (0, 0, 1)  # blue
            elif 'naive' in best['label'].lower():
                color = (1, 0, 0)  # red

            ax.plot(df_m['rc'], df_m['pr'], '--' if best['is_baseline'] else '-', label="{} (F={:.2f},C={:.2f})".format(best['label'], best['f'], best['cov_f']), color=color)
            plt.plot(best['rc'], best['pr'], 'o', color=color)

        ax.legend()
        # ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', prop={'size': 15})
        # ax.legend(bbox_to_anchor=(0, 1), loc='lower left', prop={'size': 15}, ncol=int(math.sqrt(len(labels)))-1)
        ax.set_ylabel("Precision") #, fontsize=18)
        ax.set_xlabel("Recall") #, fontsize=18)
        plt.xlim([0, 1])
        plt.ylim([0, 1])
        plt.savefig("{}/prrc_{}.png".format(out_folder, ns), bbox_inches='tight')
        plt.close('all')


def plot_mi_ru(data, out_folder):

    cmap = plt.get_cmap('tab10')
    groups = data.index.get_level_values('group').unique()

    for ns, df_ns in data.groupby(level='ns'):
        fig, ax = plt.subplots(figsize=(10, 10))

        # identify the best method for each group based on s min
        best_methods = []
        for group, df_g in df_ns[df_ns['cov_s'] > 0].groupby(level='group'):
            best_methods.append(df_g['s'].idxmin())

        # sort methods based on s
        df_best_methods = df_ns[df_ns['cov_s'] > 0].loc[best_methods].sort_values(by=['is_baseline', 's'], ascending=[True, True])

        # sort threshold
        df_best_methods = df_best_methods.set_index('tau', append=True).groupby(level='method', group_keys=False, sort=False).apply(
            lambda x: x.sort_index(level='tau', ascending=False))

        # save to file
        df_best_methods.to_csv("{}/miru_{}.tsv".format(out_folder, ns), sep='\t', float_format="%.3f")

        # plot the curves
        for method, df_m in df_best_methods.groupby(level='method', sort=False):
            best = df_m.loc[df_m['s'].idxmin()]
            color = cmap.colors[groups.get_loc(best.name[1]) % len(cmap.colors)]
            if 'blast' in best['label'].lower():
                color = (0, 0, 1)  # blue
            elif 'naive' in best['label'].lower():
                color = (1, 0, 0)  # red

            ax.plot(df_m['ru'], df_m['mi'], '--' if best['is_baseline'] else '-',
                    label="{} (S={:.2f},C={:.2f})".format(best['label'], best['s'], best['cov_s']), color=color)
            plt.plot(best['ru'], best['mi'], 'o', color=color)

        ax.legend()
        ax.set_ylabel("Misinformation")  # , fontsize=18)
        ax.set_xlabel("Remaining uncertainty")  # , fontsize=18)
        plt.xlim([0, 5])
        plt.ylim(bottom=0)
        plt.savefig("{}/miru_{}.png".format(out_folder, ns), bbox_inches='tight')
        plt.close('all')

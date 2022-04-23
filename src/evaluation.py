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


# def calc_ia(go_count_file, go_adj):
#     """
#     The negative logarithm of the conditional probability of a term being
#     annotated given that all of its parents are annotated:
#     ia(t) = -log P(t=1 | Pa(t)=1),
#     where Pa(t) is the parents of term t.
#
#     W. Clark and P. Radivojac, Information theoretic evaluation of predicted
#     ontology annotations. Bioinformatics, 2013.
#     """
#
#     go_count = parse_go_count(go_count_file)
#
#     ancestor_count = {}
#     for go_id, anchestors in go_adj.items():
#         ancestor_count.setdefault(go_id, 1)
#         for a in anchestors:
#             ancestor_count[go_id] += go_count.get(a, 0)
#
#     # Calculate information content
#     go_ia = {}
#     for go in go_adj.keys():
#         if go_count.get(go):
#             go_ia[go] = max(0.0, -math.log(float(go_count[go]) / ancestor_count[go]))
#         else:
#             go_ia[go] = 0.0
#
#     return go_ia


# def plot_prrc_curve(pr, rc, name, output_path='', index=None):
#     idx = []
#     for i in range(0,len(pr)):
#         if pr[i] > 0 and rc[i] > 0 and rc[i] < 1:
#             idx.append(i)
#     fig, ax = plt.subplots()
#     ax.plot(rc[idx], pr[idx])
#     ax.set_xlabel('Recall')
#     ax.set_ylabel('Precision')
#     ax.set_title(name)
#     plt.xlim([0,1])
#     plt.ylim([0,1])
#     if index is not None:
#         plt.plot(pr[index],rc[index], 'ro')
#     plt.savefig(output_path+name+".jpg")
#     plt.close()


# plots a graph for each different namespace containing the curves for each (or the best based on f) team's prediction
# def plot_pr_rc(pr, rc, f, cov, output_path, c_dict=None, l_dict=None, baselines=None):
#     names = list(pr.keys())
#
#     # Put the baselines at the bottom
#     if baselines:
#         names = sorted(names, key=lambda x: x in baselines, reverse=True)
#
#     namespaces = pr[names[0]].keys()
#
#     for ns in namespaces:
#         f_max = []
#         fig, ax = plt.subplots(figsize=(10, 10))
#
#         for n in names:
#             if pr[n][ns].any() and rc[n][ns].any():
#                 max_idx = np.argmax(f[n][ns])
#                 idx = []
#                 for i in range(0, len(pr[n][ns])):
#                     if pr[n][ns][i] > 0 and 0 < rc[n][ns][i] < 1:
#                         idx.append(i)
#
#                 label = l_dict.get(n, n) if l_dict is not None else n
#                 color = c_dict.get(n) if c_dict is not None else None
#
#                 if baselines and n in baselines:
#                     f_max.append(0)  # trick to put baselines at the bottom of the legend
#                     ax.plot(rc[n][ns][idx], pr[n][ns][idx], '--',
#                             label="{} F={:.2f} (C={:.2f})".format(label, f[n][ns][max_idx], cov[n][ns][max_idx]),
#                             c=color, linewidth=2)
#                 else:
#                     f_max.append(f[n][ns][max_idx])
#                     ax.plot(rc[n][ns][idx], pr[n][ns][idx],
#                             label="{} F={:.2f} (C={:.2f})".format(label, f[n][ns][max_idx], cov[n][ns][max_idx]),
#                             c=color)
#
#                 plt.yticks(np.arange(0, 1, 0.1), fontsize=16)
#                 plt.xticks(fontsize=16)
#                 plt.plot(rc[n][ns][max_idx], pr[n][ns][max_idx], 'o', c=color)
#
#
#         handles, labels = ax.get_legend_handles_labels()
#         labels2 = []
#         handles2 = []
#
#         lab_idx = list(np.argsort(f_max))
#         lab_idx.reverse()
#
#         for i in lab_idx:
#             labels2.append(labels[i])
#             handles2.append(handles[i])
#
#         ax.legend(handles2, labels2, bbox_to_anchor=(1.05, 1), loc='upper left', prop={'size': 15})
#         # Set legend on top of the plot
#         # ax.legend(handles2, labels2, bbox_to_anchor=(0, 1), loc='lower left', prop={'size': 15}, ncol=int(math.sqrt(len(labels)))-1)
#         # ax.set_title(ns, fontsize=20)
#         ax.set_ylabel("Precision", fontsize=18)
#         ax.set_xlabel("Recall", fontsize=18)
#         plt.xlim([0, 1])
#         plt.ylim([0, 1])
#         plt.savefig(output_path+"/"+ns, bbox_inches='tight')
#         plt.close('all')


# plots a graph for each different namespace containing the curves for each (or the best based on f) team's prediction
def plot_pr_rc(data, out_folder):
    # Set colors
    # c_dict = {}
    cmap = plt.get_cmap('tab20')

    # for i in range(0, len(names)):
    #     c_dict[names[i]] = cmap.colors[i % len(cmap.colors)]
    #     if 'naive' in names[i]:
    #         c_dict[names[i]] = (1, 0, 0)  # red
    #     if 'blast' in names[i]:
    #         c_dict[names[i]] = (0, 0, 1)  # blue

    groups = data.index.get_level_values('group').unique()

    # cmap.colors[i % len(cmap.colors)]

    for ns, df_ns in data.groupby(level='ns'):
        fig, ax = plt.subplots(figsize=(10, 10))

        # filter groups (1 method per group
        best_methods = []
        for group, df_g in df_ns.groupby(level='group'):
            best_methods.append(df_g['f'].idxmax())

        # plor the curves
        for group, df_g in df_ns.loc[best_methods].sort_values(by=['is_baseline', 'f'], ascending=[False, False]).groupby(level='group'):
            df_g.sort_values(by='tau', inplace=True)
            df_g.set_index('tau', append=True, inplace=True)
            best = df_g.loc[df_g['f'].idxmax()]
            print(best)

            color = cmap.colors[groups.get_loc(best.name[1]) % len(cmap.colors)]
            ax.plot(df_g['rc'], df_g['pr'], label="{} F={:.2f} (C={:.2f})".format(best['label'], best['f'], best['cov_f']), color=color)
            plt.plot(best['rc'], best['pr'], 'o', color=color)

        ax.legend()
        ax.set_ylabel("Precision", fontsize=18)
        ax.set_xlabel("Recall", fontsize=18)
        plt.xlim([0, 1])
        plt.ylim([0, 1])
        plt.savefig("{}/prrc_{}.png".format(out_folder, ns), bbox_inches='tight')
        plt.close('all')



    # for n in names:
    #     if pr[n][ns].any() and rc[n][ns].any():
    #         max_idx = np.argmax(f[n][ns])
    #         idx = []
    #         for i in range(0, len(pr[n][ns])):
    #             if pr[n][ns][i] > 0 and 0 < rc[n][ns][i] < 1:
    #                 idx.append(i)
    #
    #         label = l_dict.get(n, n) if l_dict is not None else n
    #         color = c_dict.get(n) if c_dict is not None else None
    #
    #         if baselines and n in baselines:
    #             f_max.append(0)  # trick to put baselines at the bottom of the legend
    #             ax.plot(rc[n][ns][idx], pr[n][ns][idx], '--',
    #                     label="{} F={:.2f} (C={:.2f})".format(label, f[n][ns][max_idx], cov[n][ns][max_idx]),
    #                     c=color, linewidth=2)
    #         else:
    #             f_max.append(f[n][ns][max_idx])
    #             ax.plot(rc[n][ns][idx], pr[n][ns][idx],
    #                     label="{} F={:.2f} (C={:.2f})".format(label, f[n][ns][max_idx], cov[n][ns][max_idx]),
    #                     c=color)
    #
    #         plt.yticks(np.arange(0, 1, 0.1), fontsize=16)
    #         plt.xticks(fontsize=16)
    #         plt.plot(rc[n][ns][max_idx], pr[n][ns][max_idx], 'o', c=color)
    #
    #
    #     handles, labels = ax.get_legend_handles_labels()
    #     labels2 = []
    #     handles2 = []
    #
    #     lab_idx = list(np.argsort(f_max))
    #     lab_idx.reverse()
    #
    #     for i in lab_idx:
    #         labels2.append(labels[i])
    #         handles2.append(handles[i])
    #
    #     ax.legend(handles2, labels2, bbox_to_anchor=(1.05, 1), loc='upper left', prop={'size': 15})
    #     # Set legend on top of the plot
    #     # ax.legend(handles2, labels2, bbox_to_anchor=(0, 1), loc='lower left', prop={'size': 15}, ncol=int(math.sqrt(len(labels)))-1)
    #     # ax.set_title(ns, fontsize=20)
    #     ax.set_ylabel("Precision", fontsize=18)
    #     ax.set_xlabel("Recall", fontsize=18)
    #     plt.xlim([0, 1])
    #     plt.ylim([0, 1])
    #     plt.savefig(output_path+"/"+ns, bbox_inches='tight')
    #     plt.close('all')


def plot_mi_ru():
    # _data = data[data['cov_s'] > 0]
    # best_s = _data.groupby(['ns', 'method'])['s'].transform(min) == _data['s']
    # _data[best_s].drop_duplicates().to_csv(out_folder + "/results_s.tsv",
    #                                        columns=['ns', 'method', 'group', 'label', 'mi', 'ru', 's', 'cov_s'],
    #                                        index=False, float_format='%.3f', sep='\t')

    return
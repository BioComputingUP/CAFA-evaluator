from graph import Graph, propagate
from parser import parse_obo, gt_parser, split_pred_parser, parse_ia_dict
from evaluation import get_leafs_idx, get_roots_idx, compute_metrics, compute_f, compute_s, plot_pr_rc, plot_mi_ru
import argparse
import logging
import os
import numpy as np
import pandas as pd

tau_arr = np.arange(0.01, 1, 0.01)  # Tau array, used to compute metrics at different score thresholds

parser = argparse.ArgumentParser(description='CAFA-evaluator. Calculate precision-recall plots and F-max')
parser.add_argument('obo_file', help='Ontology file in OBO format')
parser.add_argument('pred_dir', help='Predictions directory. Sub-folders are iterated recursively')
parser.add_argument('gt_file', help='Ground truth file')
parser.add_argument('-out_dir', help='Output directory. Default to \"results\" in the current directory', default='results')
parser.add_argument('-names', help='File with methods information (header: filename, group, label, is_baseline)')
parser.add_argument('-ia', help='File with information accretion (header: term, information_accretion)')
args = parser.parse_args()

obo_file = args.obo_file
pred_folder = os.path.normpath(args.pred_dir) + "/"  # add the tailing "/"
out_folder = os.path.normpath(args.out_dir) + "/"
ia_file = args.ia

# Create output folder here in order to store the log file
if not os.path.isdir(out_folder):
    os.mkdir(out_folder)

# Set the logger
logFormatter = logging.Formatter("%(asctime)s [%(levelname)-5.5s] %(message)s")
rootLogger = logging.getLogger()
rootLogger.setLevel(logging.INFO)

fileHandler = logging.FileHandler("{0}/info.log".format(out_folder))
fileHandler.setFormatter(logFormatter)
rootLogger.addHandler(fileHandler)

consoleHandler = logging.StreamHandler()
consoleHandler.setFormatter(logFormatter)
rootLogger.addHandler(consoleHandler)

# Parse the OBO file and creates a different graph for each namespace
ontologies = []
ns_dict = {}
for ns, go_dict in parse_obo(obo_file).items():
    ontologies.append(Graph(ns, go_dict))
    for term in go_dict:
        ns_dict[term] = go_dict[term]['namespace']
    logging.info("{} {} roots, {} leaves".format(ns, len(get_roots_idx(ontologies[-1].dag)), len(get_leafs_idx(ontologies[-1].dag))))


# Set prediction files
pred_files = []
for root, dirs, files in os.walk(pred_folder):
    for file in files:
        pred_files.append(os.path.join(root, file))
logging.debug("Prediction paths {}".format(pred_files))


# Set method information (optional)
methods = None
if args.names:
    methods = pd.read_csv(args.names, delim_whitespace=True)
logging.debug("Methods {}".format(methods))


# Parse and set information accretion
ia_dict = None
if ia_file is not None:
    ia_dict = parse_ia_dict(ia_file)
    for ont in ontologies:
        ont.set_ia_array(ia_dict)

# Parse ground truth file
ne = {}  # {filename: {namespace: no_proteins}}
gt = gt_parser(args.gt_file, ontologies, ns_dict)
for ns in gt:
    ont = [o for o in ontologies if o.namespace == ns][0]
    # print(gt[ns].matrix.sum(axis=0))  # statistics
    _ = propagate(gt[ns], ont, ont.order)  # propagate ancestors
    # print(gt[ns].matrix.sum(axis=0))  # statistics after ancestors
    ne[ns] = gt[ns].matrix.shape[0]  # number of proteins
    logging.info('{} proteins: {}'.format(ns, gt[ns].matrix.shape[0]))


# Parse prediction files.
# Proteins not in the ground truth are excluded as well as terms not in the ontology
# Metrics are divided only at the end of the cycle
pr_sum = {}
rc_sum = {}
wpr_sum = {}
wrc_sum = {}
mi_sum = {}
ru_sum = {}
data = {'ns': [], 'method': [], 'tau': [], 'pr': [], 'rc': [], 'f': [], 'cov_f': [], 'wpr': [], 'wrc': [], 'wf': [], 'cov_wf': [], 'mi': [], 'ru': [], 's': [], 'cov_s': []}

for file_name in pred_files:
    method = file_name.replace(pred_folder, '').replace('/', '_')
    predictions = split_pred_parser(file_name, ontologies, gt, ns_dict)

    for ns in gt:
        pr_sum[ns] = np.zeros((len(tau_arr), 2), dtype='float')
        rc_sum[ns] = np.zeros(len(tau_arr), dtype='float')

        wpr_sum[ns] = np.zeros((len(tau_arr), 2), dtype='float')
        wrc_sum[ns] = np.zeros(len(tau_arr), dtype='float')

        ru_sum[ns] = np.zeros((len(tau_arr), 2), dtype='float')
        mi_sum[ns] = np.zeros((len(tau_arr), 2), dtype='float')

    for p in predictions:
        ns = p.namespace
        if ns in gt:
            ont = [o for o in ontologies if o.namespace == ns][0]
            _ = propagate(p, ont, ont.order)

            pr, rc, wpr, wrc, ru, mi = compute_metrics(p, gt[ns], tau_arr, ont.toi, ont.ia_array)
            pr_sum[ns] += pr
            rc_sum[ns] += rc

            wpr_sum[ns] += wpr
            wrc_sum[ns] += wrc

            ru_sum[ns] += ru
            mi_sum[ns] += mi

    # compute the actual value of precision and recall for each threshold
    for ns in gt:

        data['ns'].extend([ns] * len(tau_arr))
        data['method'].extend([method] * len(tau_arr))

        data['tau'].extend(tau_arr)

        # Precision recall
        n = pr_sum[ns][:, 0]
        d = pr_sum[ns][:, 1]
        _pr = np.divide(n, d, out=np.zeros_like(n, dtype='float'), where=d != 0)

        _rc = rc_sum[ns] / ne[ns]

        data['pr'].extend(_pr)
        data['rc'].extend(_rc)
        data['f'].extend(compute_f(_pr, _rc))
        data['cov_f'].extend(d / ne[ns])

        # Mi, ru
        n = ru_sum[ns][:, 0]
        d = ru_sum[ns][:, 1]
        _ru = np.divide(n, d, out=np.zeros_like(n, dtype='float'), where=d != 0)

        n = mi_sum[ns][:, 0]
        d = mi_sum[ns][:, 1]
        _mi = np.divide(n, d, out=np.zeros_like(n, dtype='float'), where=d != 0)

        data['ru'].extend(_ru)
        data['mi'].extend(_mi)
        data['s'].extend(compute_s(_ru, _mi))
        data['cov_s'].extend(d / ne[ns])

        # Weighted precision recall
        n = wpr_sum[ns][:, 0]
        d = wpr_sum[ns][:, 1]
        _wpr = np.divide(n, d, out=np.zeros_like(n, dtype='float'), where=d != 0)

        _wrc = wrc_sum[ns] / ne[ns]

        data['wpr'].extend(_wpr)
        data['wrc'].extend(_wrc)
        data['wf'].extend(compute_f(_wpr, _wrc))
        data['cov_wf'].extend(d / ne[ns])


# Create the dataframe
data = pd.DataFrame(data)

# Add method labels and groups
if methods is None:
    group_col = 'method'
    data.set_index(['ns', 'method'], inplace=True)
else:
    group_col = 'group'
    data = pd.merge(data, methods, left_on='method', right_on='name', how='left')
    data['group'].fillna(data['method'], inplace=True)
    data['label'].fillna(data['method'], inplace=True)
    data['is_baseline'].fillna(False, inplace=True)
    data.set_index(['ns', 'group', 'method'], inplace=True)

# Groups are defined here so that they have the same colors for different namespaces
groups = data.index.get_level_values(group_col).unique()

for ns, df_ns in data.groupby(level='ns'):

    # Identify the best method for each group based on f max
    best_methods = []
    for group, df_g in df_ns.groupby(level=group_col):
        best_methods.append(df_g['f'].idxmax())

    # Sort the best methods based on f
    if methods is None:
        df_best_methods = df_ns.loc[best_methods].sort_values(by=['f'], ascending=[False])
    else:
        df_best_methods = df_ns.loc[best_methods].sort_values(by=['is_baseline', 'f'], ascending=[True, False])

    # Re-sort thresholds
    df_best_methods = df_best_methods.set_index('tau', append=True).groupby(
        level='method', group_keys=False, sort=False).apply(lambda x: x.sort_index(level='tau', ascending=False))

    # Save to file
    df_best_methods.to_csv("{}/prrc_{}.tsv".format(out_folder, ns), sep='\t', float_format="%.3f")

    # Plot
    plot_pr_rc(df_best_methods, groups, "{}/prrc_{}.png".format(out_folder, ns))

if ia_dict is not None:
    for ns, df_ns in data.groupby(level='ns'):

        # Identify the best method for each group based on f max
        best_methods = []
        for group, df_g in df_ns[df_ns['cov_s'] > 0].groupby(level='group'):
            best_methods.append(df_g['s'].idxmin())

        # Sort the best methods based on f
        if methods is None:
            df_best_methods = df_ns[df_ns['cov_s'] > 0].loc[best_methods].sort_values(by=['s'], ascending=[True])
        else:
            df_best_methods = df_ns[df_ns['cov_s'] > 0].loc[best_methods].sort_values(by=['is_baseline', 's'],
                                                                                      ascending=[True, True])
        # Re-sort thresholds
        df_best_methods = df_best_methods.set_index('tau', append=True).groupby(level='method', group_keys=False,
                                                                                sort=False).apply(
            lambda x: x.sort_index(level='tau', ascending=False))

        # Save to file
        df_best_methods.to_csv("{}/miru_{}.tsv".format(out_folder, ns), sep='\t', float_format="%.3f")

        # Plot
        plot_mi_ru(df_best_methods, groups, "{}/miru_{}.png".format(out_folder, ns))


import argparse
import logging
import os
import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt

from graph import Graph, propagate
from parser import obo_parser, gt_parser, pred_parser, ia_parser
from evaluation import get_leafs_idx, get_roots_idx, evaluate_prediction

# Tau array, used to compute metrics at different score thresholds
tau_arr = np.arange(0.01, 1, 0.01)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='CAFA-evaluator. Calculate precision-recall plots and F-max')
    parser.add_argument('obo_file', help='Ontology file in OBO format')
    parser.add_argument('pred_dir', help='Predictions directory. Sub-folders are iterated recursively. '
                                         'Files in the same sub-folder are merged')
    parser.add_argument('gt_file', help='Ground truth file')
    parser.add_argument('-out_dir', default='results',
                        help='Output directory. Default to \"results\" in the current directory')
    parser.add_argument('-ia', help='File with information accretion (header: term, information_accretion)')
    parser.add_argument('-no_orphans', action='store_true', default=False,
                        help='Consider terms without parents, e.g. the root(s)')
    parser.add_argument('-norm', choices=['cafa', 'pred', 'gt'], default='cafa',
                        help='Normalize as in CAFA. Consider predicted targets (pred). '
                             'Consider all ground truth proteins (gt)')
    parser.add_argument('-names', help='File with methods information (header: filename group label is_baseline)')
    args = parser.parse_args()

    # Create output folder here in order to store the log file
    out_folder = os.path.normpath(args.out_dir) + "/"
    if not os.path.isdir(out_folder):
        os.mkdir(out_folder)

    # Set the logger
    logFormatter = logging.Formatter("%(asctime)s [%(levelname)-5.5s] %(message)s")
    rootLogger = logging.getLogger()
    # rootLogger.setLevel(logging.DEBUG)
    rootLogger.setLevel(logging.INFO)

    fileHandler = logging.FileHandler("{0}/info.log".format(out_folder))
    fileHandler.setFormatter(logFormatter)
    rootLogger.addHandler(fileHandler)

    consoleHandler = logging.StreamHandler()
    consoleHandler.setFormatter(logFormatter)
    rootLogger.addHandler(consoleHandler)

    # Parse and set information accretion (optional)
    ia_dict = None
    if args.ia is not None:
        ia_dict = ia_parser(args.ia)

    # Parse the OBO file and creates a different graph for each namespace
    ontologies = []
    for ns, terms_dict in obo_parser(args.obo_file).items():
        ontologies.append(Graph(ns, terms_dict, ia_dict, not args.no_orphans))
        logging.info("Ontology: {}, roots {}, leaves {}".format(ns, len(get_roots_idx(ontologies[-1].dag)), len(get_leafs_idx(ontologies[-1].dag))))

    # Set prediction files
    pred_folder = os.path.normpath(args.pred_dir) + "/"  # add the tailing "/"
    pred_files = []
    for root, dirs, files in os.walk(pred_folder):
        for file in files:
            pred_files.append(os.path.join(root, file))
    logging.debug("Prediction paths {}".format(pred_files))

    # Parse ground truth file
    gt = gt_parser(args.gt_file, ontologies)

    # Parse prediction files and perform evaluation
    dfs = []
    for file_name in pred_files:
        prediction = pred_parser(file_name, ontologies, gt)
        df_pred = evaluate_prediction(prediction, gt, ontologies, tau_arr, args.norm)
        df_pred['method'] = file_name.replace(pred_folder, '').replace('/', '_')
        dfs.append(df_pred)
    df = pd.concat(dfs)

    # Add method labels and groups
    if args.names is None:
        df['group'] = df['method']
        df['label'] = df['method']
        df['is_baseline'] = False
    else:
        # Set method information (optional)
        methods = pd.read_csv(args.names, delim_whitespace=True)
        logging.debug("Methods {}".format(methods))
        df = pd.merge(df, methods, left_on='method', right_on='filename', how='left')
        df['group'].fillna(df['method'], inplace=True)
        df['label'].fillna(df['method'], inplace=True)
        if 'is_baseline' not in df:
            df['is_baseline'] = False
        else:
            df['is_baseline'].fillna(False, inplace=True)

    df = df[df['cov'] > 0].reset_index(drop=True)
    df.set_index(['group', 'label', 'ns', 'tau'], inplace=True)

    cmap = plt.get_cmap('tab20')
    df['colors'] = df.index.get_level_values('group')
    df['colors'] = pd.factorize(df['colors'])[0]
    df['colors'] = df['colors'].apply(lambda x: cmap.colors[x % len(cmap.colors)])

    for metric, cols in [('f', ['rc', 'pr']), ('wf', ['wrc', 'wpr']), ('s', ['ru', 'mi'])]:

        index_best = df.groupby(level=['group', 'ns'])[metric].idxmax() if metric in ['f', 'wf'] else df.groupby(['group', 'ns'])[metric].idxmin()
        df_best = df.loc[index_best, ['cov', 'colors'] + cols + [metric]]
        df_methods = df.reset_index('tau').loc[[ele[:-1] for ele in index_best], ['tau', 'cov', 'colors'] + cols + [metric]].sort_index()

        df_best['max_cov'] = df_methods.groupby(level=['group', 'label', 'ns'])['cov'].max()
        df_best['label'] = df_best.index.get_level_values('label')
        df_best['label'] = df_best.agg(lambda x: f"{x['label']} ({metric.upper()}={x[metric]:.3f} C={x['max_cov']:.3f})", axis=1)

        df_methods.reset_index().drop(columns=['colors']).to_csv('{}/eval_{}.tsv'.format(out_folder, metric), float_format="%.3f", sep="\t", index=False)

        df_hmean = df_best.groupby(level='group')[['cov', 'max_cov', metric]].agg(stats.hmean)
        df_hmean.to_csv('{}/hmean_{}.tsv'.format(out_folder, metric), float_format="%.3f", sep="\t")

        for ns, df_g in df_best.groupby(level='ns'):
            fig, ax = plt.subplots(figsize=(15, 15))

            for i, (index, row) in enumerate(df_g.sort_values(by=[metric, 'max_cov'], ascending=[False, False]).iterrows()):
                data = df_methods.loc[index[:-1]]
                ax.plot(data[cols[0]], data[cols[1]], color=row['colors'], label=row['label'], zorder=1000 - i)
                ax.plot(row[cols[0]], row[cols[1]], color=row['colors'], marker='o', zorder=1000 - i)

            ax.legend()

            plt.xlim(0, max(1, df_methods.loc[:, :, ns][cols[0]].max()))
            plt.ylim(0, max(1, df_methods.loc[:, :, ns][cols[1]].max()))

            plt.savefig("{}/fig_{}_{}.png".format(out_folder, metric, ns), bbox_inches='tight')
            plt.clf()

import argparse
import logging
import os
import pandas as pd
import numpy as np

from graph import Graph, propagate
from parser import obo_parser, gt_parser, pred_parser, ia_parser
from evaluation import get_leafs_idx, get_roots_idx, evaluate_prediction
from plot import get_best_methods, plot_curves

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
        df = pd.merge(df, methods, left_on='method', right_on='name', how='left')
        df['group'].fillna(df['method'], inplace=True)
        df['label'].fillna(df['method'], inplace=True)
        df['is_baseline'].fillna(False, inplace=True)

    df.set_index(['ns', 'group', 'label'], inplace=True)

    # Groups are defined here so that they will have the same color in different namespaces
    groups = df.index.get_level_values('group').unique()

    # Plot and write results to file
    for ns, df_ns in df.groupby(level='ns'):

        # Precision recall
        df_best_methods = get_best_methods(df_ns, 'f', False)
        df_best_methods.to_csv("{}/prrc_{}.tsv".format(out_folder, ns), sep='\t', float_format="%.3f",
                               columns=['pr', 'rc', 'f', 'cov'])

        plot_curves("{}/prrc_{}.png".format(out_folder, ns), df_best_methods, groups, False,
                    'f', 'rc', 'pr', 'Recall', 'Precision',
                    x_lim=(0, 1), y_lim=(0, 1),
                    label_colors={'blast': (0, 0, 1), 'naive': (1, 0, 0)})

        if args.ia is not None:

            # Weighted precision recall
            df_best_methods = get_best_methods(df_ns, 'wf', False)
            df_best_methods.to_csv("{}/wprrc_{}.tsv".format(out_folder, ns), sep='\t', float_format="%.3f",
                                   columns=['wpr', 'wrc', 'wf', 'cov'])

            plot_curves("{}/wprrc_{}.png".format(out_folder, ns), df_best_methods, groups, False,
                        'wf', 'wrc', 'wpr', 'Recall', 'Precision',
                        x_lim=(0, 1), y_lim=(0, 1),
                        label_colors={'blast': (0, 0, 1), 'naive': (1, 0, 0)})

            # Misinformation
            df_best_methods = get_best_methods(df_ns, 's', True)
            df_best_methods.to_csv("{}/miru_{}.tsv".format(out_folder, ns), sep='\t', float_format="%.3f",
                                   columns=['mi', 'ru', 's', 'cov'])

            plot_curves("{}/miru_{}.png".format(out_folder, ns), df_best_methods, groups, True,
                        's', 'ru', 'mi', 'Remaining uncertainty', 'Misinformation',
                        x_lim=(0, None), y_lim=(0, None),
                        label_colors={'blast': (0, 0, 1), 'naive': (1, 0, 0)})


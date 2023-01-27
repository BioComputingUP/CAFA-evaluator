import argparse
import logging
import os
import pandas as pd
import numpy as np
#from scipy import stats

from graph import Graph
from parser import obo_parser, gt_parser, pred_parser, ia_parser
from evaluation import get_leafs_idx, get_roots_idx, evaluate_prediction

# Tau array, used to compute metrics at different score thresholds
tau_arr = np.arange(0.01, 1, 0.01)
# tau_arr = np.arange(0.001, 1, 0.001)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='CAFA-evaluator. Calculate precision-recall curves and F-max / S-min')
    parser.add_argument('obo_file', help='Ontology file, OBO format')
    parser.add_argument('pred_dir', help='Predictions directory. Sub-folders are iterated recursively')
    parser.add_argument('gt_file', help='Ground truth file')
    parser.add_argument('-out_dir', default='results',
                        help='Output directory. By default it creates \"results/\" in the current directory')
    parser.add_argument('-ia', help='Information accretion file (columns: <term> <information_accretion>)')
    parser.add_argument('-no_orphans', action='store_true', default=False,
                        help='Consider terms without parents, e.g. the root(s), in the evaluation')
    parser.add_argument('-norm', choices=['cafa', 'pred', 'gt'], default='cafa',
                        help='Normalization strategy. i) CAFA strategy (cafa); '
                             'ii) consider predicted targets (pred); '
                             'iii) consider ground truth proteins (gt)')
    parser.add_argument('-prop', choices=['max', 'fill'], default='fill',
                        help='Ancestor propagation strategy. i) Propagate the max score of the traversed subgraph '
                             'iteratively (max); ii) Propagate with max until a different score is found (fill)')
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
        prediction = pred_parser(file_name, ontologies, gt, args.prop)
        df_pred = evaluate_prediction(prediction, gt, ontologies, tau_arr, args.norm)
        df_pred['filename'] = file_name.replace(pred_folder, '').replace('/', '_')
        dfs.append(df_pred)
        logging.info("Prediction: {}, evaluated".format(file_name))
    df = pd.concat(dfs)

    # Save the dataframe
    df = df[df['cov'] > 0].reset_index(drop=True)
    df.set_index(['filename', 'ns', 'tau'], inplace=True)
    df.to_csv('{}/evaluation_all.tsv'.format(out_folder), float_format="%.5f", sep="\t")

    # Calculate harmonic mean across namespaces for each evaluation metric
    for metric, cols in [('f', ['rc', 'pr']), ('wf', ['wrc', 'wpr']), ('s', ['ru', 'mi'])]:

        index_best = df.groupby(level=['filename', 'ns'])[metric].idxmax() if metric in ['f', 'wf'] else df.groupby(['filename', 'ns'])[metric].idxmin()

        df_best = df.loc[index_best]
        df_best['max_cov'] = df.reset_index('tau').loc[[ele[:-1] for ele in index_best]].groupby(level=['filename', 'ns'])['cov'].max()
        df_best.to_csv('{}/evaluation_best_{}.tsv'.format(out_folder, metric), float_format="%.5f", sep="\t")

        # TODO if a namespace is not predicted, create a row with all zeros
        # df_hmean = df_best.groupby(level='filename')[cols + ['cov', 'max_cov'] + [metric]].agg(stats.hmean).sort_values(metric, ascending=False if metric in ['f', 'wf'] else True)
        # df_hmean.to_csv('{}/evaluation_hmean_{}.tsv'.format(out_folder, metric), float_format="%.3f", sep="\t")


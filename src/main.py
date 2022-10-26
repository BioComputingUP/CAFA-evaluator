import argparse
import logging
import os
import pandas as pd

from graph import Graph, propagate
from parser import obo_parser, gt_parser, pred_parser, ia_parser
from evaluation import get_leafs_idx, get_roots_idx, evaluate_prediction
from plot import get_best_methods, plot_curves


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='CAFA-evaluator. Calculate precision-recall plots and F-max')
    parser.add_argument('obo_file', help='Ontology file in OBO format')
    parser.add_argument('pred_dir', help='Predictions directory. Sub-folders are iterated recursively')
    parser.add_argument('gt_file', help='Ground truth file')
    parser.add_argument('-out_dir', help='Output directory. Default to \"results\" in the current directory', default='results')
    parser.add_argument('-names', help='File with methods information (header: filename, group, label, is_baseline)')
    parser.add_argument('-ia', help='File with information accretion (header: term, information_accretion)')
    args = parser.parse_args()

    # Create output folder here in order to store the log file
    out_folder = os.path.normpath(args.out_dir) + "/"
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

    # Parse and set information accretion (optional)
    ia_dict = None
    if args.ia is not None:
        ia_dict = ia_parser(args.ia)

    # Parse the OBO file and creates a different graph for each namespace
    ontologies = []
    ns_dict = {}
    for ns, go_dict in obo_parser(args.obo_file).items():
        ontology = Graph(ns, go_dict, ia_dict)
        ontologies.append(ontology)
        for term in ontology.go_terms:
            ns_dict[term] = ontology.go_terms[term]['namespace']
        logging.info("Ontology: {}, roots {}, leaves {}".format(ns, len(get_roots_idx(ontologies[-1].dag)), len(get_leafs_idx(ontologies[-1].dag))))

    # Set prediction files
    pred_folder = os.path.normpath(args.pred_dir) + "/"  # add the tailing "/"
    pred_files = []
    for root, dirs, files in os.walk(pred_folder):
        for file in files:
            pred_files.append(os.path.join(root, file))
    logging.debug("Prediction paths {}".format(pred_files))

    # Parse ground truth file
    ne = {}  # {filename: {namespace: no_proteins}}
    gt = gt_parser(args.gt_file, ontologies, ns_dict)
    for ns in gt:
        ont = [o for o in ontologies if o.namespace == ns][0]
        # print(gt[ns].matrix.sum(axis=0))  # statistics
        _ = propagate(gt[ns], ont, ont.order)  # propagate ancestors
        # print(gt[ns].matrix.sum(axis=0))  # statistics after ancestors
        ne[ns] = gt[ns].matrix.shape[0]  # number of proteins

    # Parse prediction files and perform evaluation
    dfs = []
    for file_name in pred_files:
        prediction = pred_parser(file_name, ontologies, gt, ns_dict)
        df_pred = evaluate_prediction(prediction, gt, ontologies, ne)
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


import argparse
from cafaeval.evaluation import cafa_eval, write_results
import logging


def command_line():

    logging.info("CAFA-evaluator. Calculate precision-recall curves and F-max / S-min")

    parser = argparse.ArgumentParser(description='CAFA-evaluator. Calculate precision-recall curves and F-max / S-min')

    parser.add_argument('obo_file', help='Ontology file, OBO format')
    parser.add_argument('pred_dir', help='Predictions directory. Sub-folders are iterated recursively')
    parser.add_argument('gt_file', help='Ground truth file')

    parser.add_argument('-out_dir', default='results',
                        help='Output directory. By default it creates \"results/\" in the current directory')
    parser.add_argument('-ia', help='Information Accretion file (columns: <term> <information_accretion>)')
    parser.add_argument('-no_orphans', action='store_true', default=False,
                        help='Exclude terms without parents, e.g. the root(s), in the evaluation')
    parser.add_argument('-norm', choices=['cafa', 'pred', 'gt'], default='cafa',
                        help='cafa - implements the CAFA normalization strategy. '
                             'Precision is normalized by the number of predicted targets, '
                             'all other metrics by the number of ground truth targets. '
                             'pred - All metrics are normalized by the number of predicted targets. '
                             'gt - All metrics are normalized by the number of ground truth targets')
    parser.add_argument('-prop', choices=['max', 'fill'], default='max',
                        help='Ancestor propagation strategy. max - Propagate the max score of the traversed subgraph '
                             'iteratively. fill - Propagate with max until a different score is found')
    parser.add_argument('-th_step', type=float, default=0.01,
                        help='Threshold step size in the range [0, 1). A smaller step, means more calculation')
    parser.add_argument('-max_terms', type=int, default=None,
                        help='Number of terms for protein and namespace to consider in the evaluation')
    parser.add_argument('-threads', type=int, default=4,
                        help='Parallel threads. 0 means use all available CPU threads. '
                             'Do not use multithread if you are short in memory')
    parser.add_argument('-log_level', type=str, choices=['debug', 'info', 'warning', 'error', 'critical'],
                        default='info', help='Log level')

    args = parser.parse_args()

    # Set the logger
    logging.basicConfig()
    root_logger = logging.getLogger()  # root logger
    root_logger.setLevel(logging.getLevelName(args.log_level.upper()))
    log_formatter = logging.Formatter("%(asctime)s [%(levelname)-5.5s] %(message)s")
    root_handler = root_logger.handlers[0]
    root_handler.setFormatter(log_formatter)

    # Run the evaluation
    df, dfs_best = cafa_eval(args.obo_file, args.pred_dir, args.gt_file,
                             ia=args.ia, no_orphans=args.no_orphans, norm=args.norm, prop=args.prop,
                             max_terms=args.max_terms, th_step=args.th_step, n_cpu=args.threads)

    # Write the results
    write_results(df, dfs_best, out_dir=args.out_dir, th_step=args.th_step)


if __name__ == "__main__":

    command_line()

from graph import top_sort, propagate_gt, propagate_pred
from parser import parse_obo, parse_benchmark, gt_parser, split_pred_parser
from evaluation import get_toi_idx, get_leafs_idx, get_roots_idx, compute_f_metrics, compute_f, plot_graphs
import time
import argparse
import logging
import os
import numpy as np
import matplotlib.pyplot as plt
import sys



namespaces = {"bpo": "biological_process", "cco": "cellular_component", "mfo": "molecular_function",
 "disorderfunction": "Disorder_function", "interactionpartner": "Interaction_partner",
 "structuralstate": "Structural_state", "structuraltransition": "Structural_transition"}


parser = argparse.ArgumentParser(description='CAFA-evaluator. Calculate precision-recall plots and F-max')
parser.add_argument('obo_file', help='Ontology file in OBO format')
parser.add_argument('pred_dir', help='Predictions directory. Sub-folders are iterated recursively')
parser.add_argument('gt_dir', help='Ground truth directory. Sub-folders are iterated recursively')
parser.add_argument('-out_dir', help='Output directory. Default to \"results\"', default='results')
parser.add_argument('-filter_dir', help='Benchmark directory. By default the software consider '
                                        'all (and only) targets in the ground truth. Here you can provide'
                                        'a subset of targets to consider')
parser.add_argument('-split', help='Consider namespaces separately', default=1)
args = parser.parse_args()


tau_arr = np.arange(0.01, 1, 0.01)  # array of tau, used to compute precision and recall at different score threshold
obo_file = args.obo_file
pred_folder = os.path.normpath(args.pred_dir) + "/"  # add the tailing "/"
gt_folder = os.path.normpath(args.gt_dir) + "/"
benchmark_folder = os.path.normpath(args.filter_dir) + "/" if args.filter_dir else None
out_folder = os.path.normpath(args.out_dir) + "/"
split = args.split

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

# parsing the ontology and for each namespace in the obo file creates a different graph structure
ontologies = parse_obo(obo_file)

# Set gt files
gt_files = []
for root, dirs, files in os.walk(gt_folder):
    for file in files:
        gt_files.append(os.path.join(root, file))
logging.debug("Gt files {}".format(gt_files))

# Set prediction files
pred_files = []
for root, dirs, files in os.walk(pred_folder):
    for file in files:
        pred_files.append(os.path.join(root, file))
logging.debug("Prediction paths {}".format(pred_files))

# Set benchmark files (optional)
benchmark_files = []
if benchmark_folder:
    for root, dirs, files in os.walk(benchmark_folder):
        for file in files:
            benchmark_files.append(os.path.join(root, file))
logging.debug("Benchmark paths {}".format(benchmark_files))


pr_sum = {} 
rc_sum = {}
toi = {}
order = {}
gts = {}
toi_order = time.time()
#for each namespace initialize the precision and recall arrays and computes the terms of interest, roots and leafs 
for ont in ontologies:
    ns = ont.get_namespace()
    pr_sum[ns] = np.zeros((len(tau_arr),2), dtype='float')
    rc_sum[ns] = np.zeros(len(tau_arr), dtype='float')
    toi[ns] = get_toi_idx(ont.dag)
    roots = get_roots_idx(ont.dag)
    leafs = get_leafs_idx(ont.dag)
    """for i in roots:
        logging.info(ont.go_list[i]['id'] + " is a "+ ns +" root term") """
    """for i in leafs:
        logging.info(ont.go_list[i]['id'] + " is a "+ ns +" leaf term")"""
    order[ns] = top_sort(ont)
    logging.info("{} {} roots, {} leaves".format(ns, len(roots), len(leafs)))

# find the right ontology for each benchmark(filter) file
benchmarks = {}
gt_paths = {} 
predictions = {}
ne = {}
for k, v in namespaces.items():
    for b in benchmark_files:
        benchmarks.setdefault(v, [])
        if k.lower() in b.lower():
            benchmarks[v].append(parse_benchmark(b))
    for gt in gt_files:
        gt_paths.setdefault(v)
        if k.lower() in gt.lower():
            gt_paths[v] = gt
logging.debug("Gt paths {}".format(gt_paths))

# if a ground truth for a specific namespace is missing the namespace is ignored
unused_onts = []
for i in range(0, len(ontologies)):
    ns = ontologies[i].get_namespace()
    if gt_paths.get(ns) is None:
        unused_onts.append(i)
unused_onts.reverse()
for idx in unused_onts:
    ontologies.pop(idx)

if len(ontologies) == 0:
    raise ValueError("Ground truth filenames not matching any OBO namespace")

# parsing the ground truth for each remaining namespace
for ont in ontologies:
    ns = ont.get_namespace()
    if gt_paths.get(ns) is not None:
        b_list = benchmarks.get(ns)
        if b_list is not None:
            [b] = b_list
        else:
            b = None

        gts[ns] = gt_parser(gt_paths[ns], ont.go_terms, False, b)

        c_pre = gts[ns].matrix.sum(axis=0)

        _ = propagate_gt(gts[ns], ont, order[ns])
        if b is None:
            ne[ns] = gts[ns].matrix.shape[0]
        else:
            ne[ns] = len(b)
        logging.info("Ground truth {} targets {}".format(ns, len(gts[ns].ids)))

        c_post = gts[ns].matrix.sum(axis=0)

        for a, b, c in sorted(zip(c_post, c_pre, gts[ns].terms), reverse=True):
            print("{}\t{}\t{}".format(a - b, b, c[2]))
sys.exit(0)



# parses each prediction file, removing the proteins not contained in the ground truth
# computes precision and recall sums that will be divided only at the end of the cycle 
pr_dict = {}
rc_dict = {}
f_dict = {}
cov_dict = {}
for file_name in pred_files:
    name = file_name.replace(pred_folder, '').replace('/', '_')
    predictions = split_pred_parser(file_name, ontologies, gts)

    for p in predictions:
        ns = p.namespace
        for o in ontologies:
            if o.get_namespace() == ns:
                ont = o
        _ = propagate_pred(p, ont, order[ns])
        pr, rc = compute_f_metrics(p, gts[ns], tau_arr, toi[ns])
        pr_sum[ns] += pr
        rc_sum[ns] += rc

    # computing the actual value of precision and recall for each threshold
    pr = {}
    rc = {}
    f = {}
    cov = {}
    for ont in ontologies:
        ns = ont.get_namespace()
        n = pr_sum[ns][:, 0]
        d = pr_sum[ns][:, 1]
        pr[ns] = np.divide(n, d, out=np.zeros_like(n, dtype='float'), where=d!=0)
        rc[ns] = rc_sum[ns] / ne[ns]
        f[ns] = compute_f(pr[ns], rc[ns])
        cov[ns] = d / ne[ns]
        logging.info("F-max {} {} {}".format(ns, name, max(f[ns])))

    pr_dict[name] = pr
    rc_dict[name] = rc
    f_dict[name] = f
    cov_dict[name] = cov
    for ont in ontologies:
        ns = ont.get_namespace()
        pr_sum[ns] = np.zeros((len(tau_arr),2), dtype='float')
        rc_sum[ns] = np.zeros(len(tau_arr), dtype='float')


# file_name, group, name
groups = {
    'naive_pred_1.txt': 'naive',
    'naive_pred_2.txt': 'naive',
    'team_8_pred_1.txt': 'team_8',
    'team_4_pred_1.txt': 'team_4',
    'team_4_pred_3.txt': 'team_4',
    'team_4_pred_2.txt': 'team_4',
    'team_2_pred_1.txt': 'team_2',
    'team_2_pred_3.txt': 'team_2',
    'team_2_pred_2.txt': 'team_2',
    'team_6_pred_1.txt': 'team_6',
    'team_6_pred_2.txt': 'team_6',
    'team_1_pred_1.txt': 'team_1',
    'team_3_pred_1.txt': 'team_3',
    'team_7_pred_1.txt': 'team_7',
    'team_5_pred_1.txt': 'team_5',
    'team_5_pred_2.txt': 'team_5'
}

l_dict = {
    'naive_pred_1.txt': 'naive_m1',
    'naive_pred_2.txt': 'naive_m2',
    'team_8_pred_1.txt': 'team8_m1',
    'team_4_pred_1.txt': 'team4_m1',
    'team_4_pred_3.txt': 'team4_m3',
    'team_4_pred_2.txt': 'team4_m2',
    'team_2_pred_1.txt': 'team2_m1',
    'team_2_pred_3.txt': 'team2_m3',
    'team_2_pred_2.txt': 'team2_m2',
    'team_6_pred_1.txt': 'team6_m1',
    'team_6_pred_2.txt': 'team6_m2',
    'team_1_pred_1.txt': 'team1_m1',
    'team_3_pred_1.txt': 'team3_m1',
    'team_7_pred_1.txt': 'team7_m1',
    'team_5_pred_1.txt': 'team5_m1',
    'team_5_pred_2.txt': 'team5_m2'
}

# Write to file
names = list(pr_dict.keys())
namespaces = pr_dict[names[0]].keys()
bests = {}
with open(out_folder + "/results.tsv", "wt") as out_file:
    out_file.write("\t".join(["method", "namespace", "fmax", "precision", "recall", "coverage", "tau"]) + "\n")
    for ns in namespaces:
        data = {}
        for n in names:
            if pr_dict[n][ns].any() and rc_dict[n][ns].any():
                max_idx = np.argmax(f_dict[n][ns])
                group = groups.get(n, n)
                data.setdefault(group, []).append([n, ns, f_dict[n][ns][max_idx], pr_dict[n][ns][max_idx], rc_dict[n][ns][max_idx], cov_dict[n][ns][max_idx], tau_arr[max_idx]])

        # filter by group
        for g in data:
            _n, _ns, _fmax, _pr, _rc, _cov, _tau = sorted(data[g], key=lambda x: x[2], reverse=True)[0]
            # write to file
            out_file.write("{}\t{}\t{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}\t{:.2f}\n".format(l_dict.get(_n, _n), _ns, _fmax, _pr, _rc, _cov, _tau))

            bests.setdefault(_n, {}).setdefault(_ns, None)

# Remove non-best methods
for n in list(pr_dict.keys()):
    if n not in bests:
        del pr_dict[n]
        del rc_dict[n]
        del f_dict[n]
        del cov_dict[n]
    else:
        for ns in list(pr_dict[n].keys()):
            if ns not in bests[n]:
                del pr_dict[n][ns]
                del rc_dict[n][ns]
                del f_dict[n][ns]
                del cov_dict[n][ns]
        if not pr_dict[n]:
            del pr_dict[n]
            del rc_dict[n]
            del f_dict[n]
            del cov_dict[n]

# Set colors
c_dict = {}
cmap = plt.get_cmap('tab20')

for i in range(0, len(names)):
    c_dict[names[i]] = cmap.colors[i % len(cmap.colors)]
    if 'naive' in names[i]:
        c_dict[names[i]] = (0, 0, 0)

plot_graphs(pr_dict, rc_dict, f_dict, cov_dict, out_folder, c_dict, l_dict)





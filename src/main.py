from graph import top_sort, propagate_gt, propagate_pred
from parser import parse_obo, parse_benchmark, gt_parser, split_pred_parser
from evaluation import get_toi_idx, get_leafs_idx, get_roots_idx, compute_f_metrics, compute_f, plot_graphs
import time
import argparse
import logging
import os
import numpy as np


#begin = time.time()
namespaces = {"bpo": "biological_process", "cco": "cellular_component", "mfo": "molecular_function",
 "disorderfunction": "Disorder function", "interactionpartner": "Interaction partner",
 "structuralstate": "Structural state", "structuraltransition": "Structural transition"}

parser = argparse.ArgumentParser(description='CAFA-evaluator. Calculate precision-recall plots and F-max')
parser.add_argument('obo_file', help='Ontology file in OBO format')
parser.add_argument('pred_folder', help='Predictions folder. Sub-folders are iterated recursively')
parser.add_argument('gt_folder', help='Ground truth folder. Sub-folders are iterated recursively')
parser.add_argument('-output_path', help='Output directory', default='eval_results')
parser.add_argument('-split', help='Consider namespaces separately', default=1)
parser.add_argument('-gt_filters', help='Benchmark folder. By default the software consider '
                                        'all (and only) targets in the ground truth. Here you can provide'
                                        'a subset of targets to consider')
args = parser.parse_args()


tau_arr = np.arange(0.01, 1, 0.01) # array of tau, used to compute precision and recall at different score threshold
obo_file = args.obo_file
pred_folder = os.path.normpath(args.pred_folder) + "/"  # add the tailing "/"
gt_folder = os.path.normpath(args.gt_folder) + "/"
split = args.split
benchmark_folder = args.gt_filters
output_path = args.output_path 

if not os.path.isdir(output_path):
    os.mkdir(output_path)

# Set the logger
logFormatter = logging.Formatter("%(asctime)s [%(levelname)-5.5s] %(message)s")
rootLogger = logging.getLogger()
rootLogger.setLevel(logging.INFO)

fileHandler = logging.FileHandler("{0}/info.log".format(output_path))
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
logging.info("Prediction paths {}".format(pred_files))

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
#print("toi and order:", time.time()-toi_order)

# find the right ontology for each benchmark(filter) file
benchmarks = {}
gt_paths = {} 
predictions = {}
ne = {}
# gt_pat_set = time.time()
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
#print("setting gt path:", time.time()-gt_pat_set)

# if a ground truth for a specific namespace is missing the namespace is ignored
# rm_onts = time.time()
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

#print("removing unused ontologies:", time.time()-rm_onts)


# parsing the ground truth for each remaining namespace
# gt_parsing = time.time()
for ont in ontologies:
    ns = ont.get_namespace()
    if gt_paths.get(ns) is not None:
        b_list = benchmarks.get(ns)
        if b_list is not None:
            [b] = b_list
        else:
            b = None

        gts[ns] = gt_parser(gt_paths[ns], ont.go_terms, False, b)
        _ = propagate_gt(gts[ns], ont, order[ns])
        if b is None:
            ne[ns] = gts[ns].matrix.shape[0]
        else:
            ne[ns] = len(b)
        logging.info("Ground truth {} targets {}".format(ns, len(gts[ns].ids)))
#print("gt parsing and propagation:", time.time()-gt_parsing)

# parses each prediction file, removing the proteins not contained in the ground truth
# computes precision and recall sums that will be divided only at the end of the cycle 
pr_dict = {}
rc_dict = {}
f_dict = {}
# pred_and_f = time.time()
for file_name in pred_files:
    name = file_name.replace(pred_folder, '').replace('/', '_')
    # pred_parse_create = time.time()
    predictions = split_pred_parser(file_name, ontologies, gts)
    #print("prediction parse and creation:", time.time()-pred_parse_create)
    for p in predictions:
        ns = p.namespace
        for o in ontologies:
            if o.get_namespace() == ns:
                ont = o
        # propagation = time.time()
        _ = propagate_pred(p, ont, order[ns])
        #print("propagation:", time.time()-propagation)
        # computing_prrc = time.time()
        pr, rc = compute_f_metrics(p, gts[ns], tau_arr, toi[ns])
        pr_sum[ns] += pr
        rc_sum[ns] += rc
        #print("computing pr and rc terms:", time.time()-computing_prrc, ns)

    # computing the actual value of precision and recall for each threshold
    pr = {}
    rc = {}
    f = {}
    for ont in ontologies:
        ns = ont.get_namespace()
        n = pr_sum[ns][:,0]
        d = pr_sum[ns][:,1]
        pr[ns] = np.divide(n, d, out=np.zeros_like(n, dtype='float'), where=d!=0)
        rc[ns] = rc_sum[ns] / ne[ns]
        f[ns] = compute_f(pr[ns], rc[ns])
        logging.info("F-max {} {} {}".format(ns, name, max(f[ns])))

    pr_dict[name] = pr
    rc_dict[name] = rc
    f_dict[name] = f
    for ont in ontologies:
        ns = ont.get_namespace()
        pr_sum[ns] = np.zeros((len(tau_arr),2), dtype='float')
        rc_sum[ns] = np.zeros(len(tau_arr), dtype='float')
    #print("time for single prediction:", time.time()-for_time)
#print("parsing predictions and computing f:", time.time()-pred_and_f)

plot_graphs(pr_dict, rc_dict, f_dict, output_path, tau_arr)

#print("Time elapsed:", time.time()-begin)






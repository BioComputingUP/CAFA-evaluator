from graph import *
from parser import * 
from evaluation import *
import time
import argparse
import logging
from os import listdir
import os

#take a prediction file containing paths to prediction file cycling over it calculate pr and rc 
#take a gt file, this will be reused for each prediction



# input FILE_PREDIZIONE, FILE_GT, FILE_OBO, array di tau, split o no, benchmark (partial evaluation) 
# output: fmax , grafico roc  
begin = time.time()
namespaces = {"bpo": "biological_process", "cco": "cellular_component", "mfo": "molecular_function",
 "disorderfunction": "Disorder function", "interactionpartner": "Interaction partner",
 "structuralstate": "Structural state", "structuraltransition": "Structural transition"}

parser = argparse.ArgumentParser(description='CAFA evaluator')
parser.add_argument('obo_file', help='File containing the ontology')
parser.add_argument('pred_folder', help='Folder containing the prediction files') 
parser.add_argument('gt_folder', help='Folder containing the ground truth files')
parser.add_argument('-output_path', help='Output directory', default='eval_results')
parser.add_argument('-split', help='0 if not splitting, 1 if splitting the ontology (put it to 1 if the ontology has separated sub-graphs', default=1)
parser.add_argument('-gt_filters', help='Folder containing the path of benchmark files to use during the evaluation') 
args = parser.parse_args()


tau_arr = np.arange(0.01, 1, 0.01) # array of tau, used to compute precision and recall at different score threshold
obo_file = args.obo_file
pred_folder = args.pred_folder
gt_folder = args.gt_folder
split = args.split
benchmark_folder = args.gt_filters
output_path = args.output_path 

if not os.path.isdir(output_path):
    os.mkdir(output_path)

logging.basicConfig(filename=output_path+"/info.log", filemode='w', format='%(message)s', level=logging.INFO)

preds = {}
for root, dirs, files in os.walk(pred_folder):
    if len(files) > 0:
        preds[root] = files

if benchmark_folder:
    benchmark_paths = [benchmark_folder+'/'+f for f in listdir(benchmark_folder)]
else:
    benchmark_paths = []
gt_files = [gt_folder+"/"+f for f in listdir(gt_folder)]


# parsing the ontology and for each namespace in the obo file creates a different graph structure
ont_time = time.time()
ontologies = parse_obo(obo_file)
#print("ontology parse:", time.time()-ont_time)
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
    logging.info("in the "+ ns +" ontology there is a total of " + str(len(roots)) + " roots")
    """for i in leafs:
        logging.info(ont.go_list[i]['id'] + " is a "+ ns +" leaf term")"""
    logging.info("in the "+ ns +" ontology there is a total of " + str(len(leafs)) + " leafs")
    order[ns] = top_sort(ont)
#print("toi and order:", time.time()-toi_order)

# find the right onotlogy for each benchmark(filter) file 
benchmarks = {}
gt_paths = {} 
predictions = {}
ne = {}
gt_pat_set = time.time()
for k, v in namespaces.items():
    for b in benchmark_paths:
        benchmarks.setdefault(v,[])
        if k.lower() in b.lower():
            benchmarks[v].append(parse_benchmark(b))
    for gt in gt_files:
        gt_paths.setdefault(v, None)
        if k.lower() in gt.lower():
            gt_paths[v] = gt 
#print("setting gt path:", time.time()-gt_pat_set)

# if a ground truth for a specific namespace is missing the namespace is ignored
rm_onts = time.time()
unused_onts = []
for i in range(0,len(ontologies)):
    ns = ontologies[i].get_namespace()
    if gt_paths.get(ns) is None:
        unused_onts.append(i)
unused_onts.reverse()
for idx in unused_onts:
    ontologies.pop(idx)
#print("removing unused ontologies:", time.time()-rm_onts)

# parsing the ground truth for each remaining namespace
gt_parsing = time.time()
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
#print("gt parsing and propagation:", time.time()-gt_parsing)

# parses each prediction file, removing the proteins not contained in the groud truth
# computes precision and recall sums that will be divided only at the end of the cycle 
pr_dict = {}
rc_dict = {}
f_dict = {}   
pred_and_f = time.time()    
for k, pred_paths in preds.items():
    for_time = time.time()
    name = k.replace(pred_folder+'/', '').replace('/','_')
    for path in pred_paths:
        #print(path)
        pred_parse_create = time.time()
        predictions = split_pred_parser(k+'/'+path, ontologies, gts)
        #print("prediction parse and creation:", time.time()-pred_parse_create)
        for p in predictions:
            ns = p.namespace
            for o in ontologies:
                if o.get_namespace() == ns:
                    ont = o
            propagation = time.time()
            _ = propagate_pred(p, ont, order[ns])
            #print("propagation:", time.time()-propagation)
            computing_prrc = time.time()
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
        print("F max "+ns.replace('_',' ')+' '+name+":", max(f[ns]))

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

            




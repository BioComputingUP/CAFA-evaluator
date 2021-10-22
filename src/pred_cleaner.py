from parser import *
import os
import argparse
from os import listdir
from graph import *
from evaluation import *

namespaces = {"bpo": "biological_process", "cco": "cellular_component", "mfo": "molecular_function",
 "disorderfunction": "Disorder function", "interactionpartner": "Interaction partner",
 "structuralstate": "Structural state", "structuraltransition": "Structural transition"}


def clean_pred(pred_file, gts, output_path):
    if not os.path.isdir(output_path):
        os.mkdir(output_path)
    rel_list = []
    namespaces = []
    for ont in ontologies:
        ns = ont.ontology_type
        namespaces.append(ns)
    parsing = time.time()
    with open(pred_file) as f:
        for line in f:
            line = line.strip().split()
            if line and len(line) > 2:
                p_id, go_id, prob = line[:3]
                if p_id.startswith("T"):
                    for ns in namespaces:
                        if p_id in gts[ns].ids:
                            rel_list.append([p_id, go_id, prob])
                            break

    filename = pred_file.split('.')[0].split('/')[-1]
    f = open(output_path + '/' +filename, 'w')
    for p_id, go_id, score in rel_list:
        f.write(p_id + ' ' + go_id + ' ' + score + '\n')

parser = argparse.ArgumentParser(description='CAFA evaluator')
parser.add_argument('obo_file', help='File containing the ontology')
parser.add_argument('pred_folder', help='Folder containing the prediction files') # pred dir oppure singolo file 
parser.add_argument('gt_folder', help='Folder containing the ground truth files')
parser.add_argument('output_path', help='Output directory')
parser.add_argument('-split', help='0 if not splitting, 1 if splitting the ontology (put it to 1 if the ontology has separated sub-graphs', default=0)
#parser.add_argument('-full_evaluation', help='1 for full evaluation, 0 for partial evaluation', default=0)
parser.add_argument('-gt_filters', help='Folder containing the path of benchmark files to use during the evaluation') #facoltativo e non interessante in realtà
args = parser.parse_args()

obo_file = args.obo_file
pred_folder = args.pred_folder
gt_folder = args.gt_folder
split = args.split
benchmark_folder = args.gt_filters
output_path = args.output_path

preds = {}
for root, dirs, files in os.walk(pred_folder):
    if len(files) > 0:
        preds[root] = files

if benchmark_folder:
    benchmark_paths = [benchmark_folder+'/'+f for f in listdir(benchmark_folder)]
else:
    benchmark_paths = []
gt_files = [gt_folder+"/"+f for f in listdir(gt_folder)]

toi = {}
order = {}
gts = {}
ontologies = parse_obo_smart(obo_file)
toi_order = time.time()
for ont in ontologies:
    ns = ont.get_namespace()
    toi[ns] = get_root_term_idx(ont.dag)
    order[ns] = top_sort(ont)

benchmarks = {}
gt_paths = {} 
predictions = {}
ne = {}
gt_pat_set = time.time()
for k, v in namespaces.items():
    for b in benchmark_paths:
        benchmarks.setdefault(v,[])
        if k in b:
            benchmarks[v].append(parse_benchmark(b))
    for gt in gt_files:
        gt_paths.setdefault(v, None)
        if k in gt:
            gt_paths[v] = gt 

unused_onts = []
for i in range(0,len(ontologies)):
    ns = ontologies[i].get_namespace()
    if gt_paths.get(ns) is None:
        unused_onts.append(i)
unused_onts.reverse()
for idx in unused_onts:
    ontologies.pop(idx)

for ont in ontologies:
    ns = ont.get_namespace()
    if gt_paths.get(ns) is not None:
        b_list = benchmarks.get(ns)
        if b_list is not None:
            [b] = b_list
        else:
            b = None
    
        gts[ns] = gt_parser2(gt_paths[ns], ont.go_terms, False, b)
        _ = propagate_gt(gts[ns], ont, order[ns])
        if b is None:
            ne[ns] = gts[ns].matrix.shape[0]
        else:
            ne[ns] = len(b)

for k, pred_paths in preds.items():
    for_time = time.time()
    name = k.replace(pred_folder+'/', '').replace('/','_')
    for path in pred_paths:
        pred_parse_create = time.time()
        predictions = clean_pred(k+'/'+path, gts, output_path)



        
    

    
    
    
from graph import *
import copy
import logging


# parses a obo file and returns a list of ontologies, one for each different namespace
def parse_obo(obo_file, rel=["is_a", "part_of"]):
    namespaces = set()
    rel_list = []
    go_dict = {}
    temp_id = ''
    temp_namespace = '' 
    temp_name = ''
    temp_def = ''
    temp_alt_id = []
    temp_rel = []
    with open(obo_file) as f:
        for line in f:
            line = line.strip().split(": ")
            if line and len(line) > 1:
                k, v = line[:2]
                if k == "id":
                    # when a new id is found first we have to input the entry in the go term dict
                    # also we will have to input a new entry when we reach EOF 
                    if temp_id != '':
                        go_dict[temp_id] = {'name': temp_name, 'namespace': temp_namespace, 'def': temp_def, 'alt_id': temp_alt_id, 'rel': temp_rel}
                    temp_alt_id = []
                    temp_id = v
                    temp_rel = []
                elif k == "alt_id":
                    temp_alt_id.append(v)
                elif k == "name":
                    temp_name = v
                elif k == "namespace" and v != 'external':
                    temp_namespace = v
                    namespaces.add(temp_namespace)
                elif k == "def":
                    temp_def = v
                elif k == "is_a" and k in rel:
                    # add (temp_id,v) tuple to the relation list
                    s = v.split('!')[0].strip()
                    # rel_list.append([temp_id, s, temp_namespace])

                    temp_rel.append(s)
                elif k == "relationship" and v.startswith("part_of") and "part_of" in rel:
                    s = v.split()[1].strip()
                    temp_rel.append(s)
                    # rel_list.append([temp_id, s, temp_namespace])

        # Last record
        go_dict[temp_id] = {'name': temp_name, 'namespace': temp_namespace, 'def': temp_def, 'alt_id': temp_alt_id, 'rel': temp_rel}

    go_struct = {}  # {ns: {term: {index: , name: , namespace: , def: }}
    go_lists = {}  # {ns: [{id: term, name:, namespace: , def:, adg: [], children: []}, ...]}
    dags = {}  # {ns: [[], ...]}  terms X terms matrix
    idxs = {}  # {ns: <counter>}
    
    for go_id, term in go_dict.items():
        ns = term['namespace']

        rel_list.extend([[go_id, rel, term['namespace']] for rel in term['rel']])

        go_lists.setdefault(ns, [])
        go_lists[ns].append({'id': go_id, 'name': term['name'], 'namespace': ns, 'def': term['def'], 'adj': [], 'children': []})

        idxs.setdefault(ns, 0)
        go_struct.setdefault(ns, {})
        go_struct[ns][go_id] = {'index': idxs[ns], 'name': term['name'], 'namespace': ns, 'def': term['def']}
        for a_id in term['alt_id']:
            go_struct.setdefault(ns, {})
            go_struct[ns][a_id] = copy.copy(go_struct[ns][go_id]) # {'index': idxs[ns], 'name': name, 'namespace': ns, 'def': definition}

        idxs[ns] += 1

    for ns in namespaces:
        go_struct[ns]['n'] = idxs[ns]
        dags[ns] = np.zeros((idxs[ns], idxs[ns]), dtype='bool')
    
    for id1, id2, ns in rel_list:
        i = go_struct[ns][id1]['index']
        j = go_struct[ns][id2]['index']
        dags[ns][i, j] = 1
        go_lists[ns][i]['adj'].append(j)
        go_lists[ns][j]['children'].append(i)

    ontologies = []
    for ns in namespaces:
        ontologies.append(Graph(dags[ns], go_struct[ns], ns, go_lists[ns]))
    
    return ontologies


# parses a file containing predictions and returns a prediction object
def pred_parser(pred_file, go_terms, gt_ids):
    idx = 0
    pr_id = {}
    rel_list = []
    with open(pred_file) as f:
        for line in f:
            line = line.strip().split()
            if line and len(line) > 2:
                p_id, go_id, prob = line[:3] 
                if p_id.startswith("T") and go_id.startswith("GO:"):
                    if pr_id.get(p_id) is None and gt_ids.get(p_id) is not None:
                        pr_id[p_id] = idx
                        idx += 1
                    rel_list.append([p_id, go_id, prob])
    n_pr = idx
    n_go = go_terms['n']
    scores = np.zeros((n_pr, n_go), dtype='float')
    for p_id, go_id, prob in rel_list:
        g = go_terms.get(go_id)
        p = pr_id.get(p_id)
        if g is not None and p is not None:
            i = p
            j = g['index']
            scores[i, j] = prob
    return Prediction(pr_id, scores, idx)


# Parses a prediction file and returns a Prediction object for each ontology in ontologies
def split_pred_parser(pred_file, ontologies , gts):
    rel_list = []
    unused_pid = set()
    tot_ids = {}
    idx = {}
    ids = {}
    namespaces = []
    for ont in ontologies:
        ns = ont.ontology_type
        idx.setdefault(ns, 0)
        ids.setdefault(ns, {})
        namespaces.append(ns)
        tot_ids.update(gts[ns].ids)

    with open(pred_file) as f:
        for line in f:
            line = line.strip().split()
            if line and len(line) > 2:
                p_id, go_id, prob = line[:3]
                if p_id.startswith("T"):
                    for ns in namespaces: 
                        if p_id not in ids[ns] and p_id in gts[ns].ids:
                            ids[ns][p_id] = idx[ns]
                            idx[ns] += 1 
                            break
                    if p_id in tot_ids:
                        rel_list.append([p_id, go_id, prob])
                    else:
                        unused_pid.add(p_id)

    # for p_id in unused_pid:
    #     logging.debug("Not in ground truth {}".format(p_id))

    logging.debug("Skipped targets: {} (not gt) {} (not pred) ({})".format(len(unused_pid),
                                                                          len(set(gts[ns].ids) - set(ids[ns].keys())),
                                                                          pred_file))
    for ns in gts:
        for p_id in set(gts[ns].ids) - set(ids[ns].keys()):
            logging.debug("Not predicted {} {}".format(ns, p_id))
            # print("Not predicted {} {}".format(ns, p_id))

    scores = {}
    np_ = {}
    nt = {}
    for ont in ontologies:
        ns = ont.ontology_type
        np_[ns] = idx[ns]
        nt[ns] = ont.go_terms['n']
        scores[ns] = np.zeros((np_[ns], nt[ns]), dtype='float')

    for p_id, go_id, prob in rel_list:
        for ont in ontologies:
            ns = ont.ontology_type
            if p_id in ids[ns] and go_id in ont.go_terms:
                i = ids[ns].get(p_id)
                j = ont.go_terms.get(go_id)['index']
                scores[ns][i, j] = prob
                break

    predictions = []
    for ns in namespaces:
        predictions.append(Prediction(ids[ns], scores[ns], idx[ns], ns))
        logging.debug("Used targets {} {} ({})".format(ns, len(ids[ns]), pred_file))
    return predictions


# parses a ground truth file, using the dictionary of go terms to remove go terms that are not contained in the ontology
# it is possible to use a benchmark file in order to filter the protein ids in the ground truth
def gt_parser(gt_file, go_terms):
    idx = 0
    rel_list = []
    gt_ids = {}
    with open(gt_file) as f:
        for line in f:
            line = line.strip().split()
            if line and len(line) == 2:
                p_id, go_id = line[:2]
                #if p_id.startswith("T"):
                if p_id not in gt_ids:
                    gt_ids[p_id] = idx
                    idx += 1
                rel_list.append([p_id, go_id])
    n = len(gt_ids)

    terms = sorted([(v['index'], k, v['name']) for k, v in go_terms.items() if k != 'n'])
    gt = np.zeros((n, go_terms['n']), dtype='bool')
    for p_id, go_id in rel_list:
        if go_id in go_terms and p_id in gt_ids:
            gt[gt_ids[p_id], go_terms[go_id]['index']] = 1
        elif go_id not in go_terms:
            logging.debug("the term " + go_id + " is not contained in the ontology")
    return GroundTruth(gt_ids, gt, terms=terms)
    

def parse_pred_paths(paths_file):
    paths = []
    with open(paths_file) as f:
        for line in f:
            paths.append(line.strip('\n'))
    return paths


# parses a benchmark file
def parse_benchmark(file):
    p_ids = {}
    with open(file) as f:
        for line in f:
            p_ids.setdefault(line.strip('\n'), 1)
    return p_ids


def parse_ia_dict(file):
    ia_dict = {}
    with open(file) as f:
        for line in f:
            if line:
                term, ia = line.strip().split()
                ia_dict[term] = float(ia)
    return ia_dict


# TOY TEST
"""[ont] = parse_obo("src/go.obo")
pred = pred_parser("src/pred.txt", ont.go_terms)
gt = gt_parser("src/gt.txt", pred, ont.go_terms)
print(ont.go_list)


print("DAG\n", ont.dag)
print("pred\n", pred.scores)
print("gt\n", gt)

print("propagated pred: \n", propagate_pred(pred.scores, ont))
print("propagated gt;\n", propagate_pred(gt, ont))   """


def mistake_counter(a1, a2):
    #print(len(a1), len(a2))
    print(a1.sum(), a2.sum())
    n = len(a1)
    indexes = []
    mis = np.zeros(n, dtype='bool')
    for i in range(0,len(a1)):
        if a1[i] != a2[i]: #a1[i] == 0 and a2[i] == 1:
            mis[i] = 1
            indexes.append(i)
    return mis.sum(), indexes

"""# TEST
[ont, _, _] = parse_obo("src/go_term.obo", True)
print(1)
pred = pred_parser("src/pred.txt", ont.go_terms)
pred2 = pred_parser("src/pred.txt", ont.go_terms)
print(2)
leaf_gt = gt_parser("src/leafonly_BPO.txt", pred, ont.go_terms)
print(leaf_gt.shape)
prop_gt = gt_parser("src/propagated_BPO.txt", pred2, ont.go_terms)
print(prop_gt.shape)
start = time.time()
gt, order = propagate_pred(leaf_gt, ont)
end = time.time()
print("Shape:",gt.shape)
print("Any:", gt.any())
print("Equal:", np.array_equal(gt[4:5,:],prop_gt[4:5,:]))
#rows1 , cols1 = np.where(gt!=0)
#rows2, cols2 = np.where(prop_gt!=0)
#print(len(rows1),len(rows2), len(cols1), len(cols2))
m, idxs = mistake_counter(gt[3,:], prop_gt[3,:])
print("mistakes:", m)
print("order:", order)
print("len order", len(order))
ids = []
for i in idxs:
    ids.append(ont.go_lis[i]['id'])
print(ids)
print(ont.dag[ont.go_terms['GO:0009991']['index'],ont.go_terms['GO:0009605']['index']])
print(pred.ids['T100900000045'])
print("time:", end-start)
#print(gt[:10,:10],prop_gt[:10,:10]) """                  
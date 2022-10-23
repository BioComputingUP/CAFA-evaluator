from graph import Prediction, GroundTruth
import numpy as np
import logging
import xml.etree.ElementTree as ET


# TODO replace the first part with a OBO parser library and ove the rest inside the Graph class
# parses a obo file and returns a list of ontologies, one for each different namespace
def parse_obo(obo_file, rel=["is_a", "part_of"]):
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
                        go_dict.setdefault(temp_namespace, {})[temp_id] = {'name': temp_name, 'namespace': temp_namespace, 'def': temp_def, 'alt_id': temp_alt_id, 'rel': temp_rel}
                    temp_alt_id = []
                    temp_id = v
                    temp_rel = []
                elif k == "alt_id":
                    temp_alt_id.append(v)
                elif k == "name":
                    temp_name = v
                elif k == "namespace" and v != 'external':
                    temp_namespace = v
                elif k == "def":
                    temp_def = v
                elif k == "is_a" and k in rel:
                    s = v.split('!')[0].strip()
                    temp_rel.append(s)
                elif k == "relationship" and v.startswith("part_of") and "part_of" in rel:
                    s = v.split()[1].strip()
                    temp_rel.append(s)

        # Last record
        go_dict.setdefault(temp_namespace, {})[temp_id] = {'name': temp_name, 'namespace': temp_namespace, 'def': temp_def, 'alt_id': temp_alt_id, 'rel': temp_rel}
    return go_dict


def gt_parser(gt_file, ontologies, ns_dict):
    """
    Parse ground truth file
    Discard terms not included in the ontology
    """
    c = 0
    gt_dict = {}
    with open(gt_file) as f:
        for line in f:
            line = line.strip().split()
            if line:
                p_id, go_id = line[:2]
                if go_id in ns_dict:
                    gt_dict.setdefault(ns_dict[go_id], {}).setdefault(p_id, []).append(go_id)
                else:
                    logging.debug("Term {} not in current ontology".format(go_id))
                c += 1
                if c == 1000:
                    break
    gts = {}
    for ont in ontologies:
        terms = sorted([(v['index'], k, v['name']) for k, v in ont.go_terms.items()])
        matrix = np.zeros((len(gt_dict[ont.namespace]), ont.idxs), dtype='bool')
        ids = {}
        for i, p_id in enumerate(gt_dict[ont.namespace]):
            ids[p_id] = i
            for go_id in gt_dict[ont.namespace][p_id]:
                matrix[i, ont.go_terms[go_id]['index']] = 1
        gts[ont.namespace] = GroundTruth(ids, matrix, terms=terms)

    return gts


# TODO improve reporting of missing targets, etc.
# Parses a prediction file and returns a Prediction object for each ontology in ontologies
def split_pred_parser(pred_file, ontologies, gts, ns_dict):

    pred_dict = {}
    with open(pred_file) as f:
        for line in f:
            line = line.strip().split()
            if line and len(line) > 2:
                p_id, go_id, prob = line[:3]
                if go_id in ns_dict and p_id in gts[ns_dict[go_id]].ids:
                    pred_dict.setdefault(ns_dict[go_id], {}).setdefault(p_id, []).append((go_id, prob))

    predictions = []
    for ont in ontologies:
        ns = ont.namespace
        if ns in pred_dict:
            matrix = np.zeros((len(pred_dict[ns]), ont.idxs), dtype='float')
            ids = {}
            for i, p_id in enumerate(pred_dict[ns].keys()):
                ids[p_id] = i
                for go_id, prob in pred_dict[ns][p_id]:
                    j = ont.go_terms.get(go_id)['index']
                    matrix[i, j] = prob
            predictions.append(Prediction(ids, matrix, len(pred_dict[ns]), ns))

    return predictions


def parse_ia_dict(file):
    ia_dict = {}
    with open(file) as f:
        for line in f:
            if line:
                term, ia = line.strip().split()
                ia_dict[term] = float(ia)
    return ia_dict


def parse_sprot(input_file, output_file):
    """
    Parse the Swiss-Prot XML annotation file
    and write a TSV file with:
    - accession
    - GO term
    - evidence code (ECO)
    """
    namespaces = {'uniprot': 'http://uniprot.org/uniprot'}
    nsl = len(namespaces['uniprot']) + 2
    with open(input_file) as f:
        with open(output_file, 'w') as fout:
            for event, elem in ET.iterparse(f, events=('start', 'end')):
                if event == 'end':
                    if elem.tag[nsl:] == 'entry':
                        acc = elem.find('uniprot:accession', namespaces).text
                        for el in elem.iterfind('uniprot:dbReference', namespaces):
                            if el.attrib['type'] == 'GO':
                                for at in el.iter():
                                    if at.attrib['type'] == 'evidence':
                                        fout.write('{}\t{}\t{}\n'.format(acc, el.attrib['id'], at.attrib['value']))
                        elem.clear()
    return

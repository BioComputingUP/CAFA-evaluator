from graph import Prediction, GroundTruth
import numpy as np
import logging
import xml.etree.ElementTree as ET


def obo_parser(obo_file, valid_rel=("is_a", "part_of")):
    """
    Parse a OBO file and returns a list of ontologies, one for each namespace.
    Obsolete terms are excluded as well as external namespaces
    """
    go_dict = {}
    term_id = None
    namespace = None
    name = None
    term_def = None
    alt_id = []
    rel = []
    obsolete = True
    with open(obo_file) as f:
        for line in f:
            line = line.strip().split(": ")
            if line and len(line) > 1:
                k = line[0]
                v = ": ".join(line[1:])
                if k == "id":
                    # Populate the dictionary with the previous entry
                    if term_id is not None and obsolete is False and namespace is not None:
                        go_dict.setdefault(namespace, {})[term_id] = {'name': name,
                                                                       'namespace': namespace,
                                                                       'def': term_def,
                                                                       'alt_id': alt_id,
                                                                       'rel': rel}
                    # Assign current term ID
                    term_id = v

                    # Reset optional fields
                    alt_id = []
                    rel = []
                    obsolete = False
                    namespace = None

                elif k == "alt_id":
                    alt_id.append(v)
                elif k == "name":
                    name = v
                elif k == "namespace" and v != 'external':
                    namespace = v
                elif k == "def":
                    term_def = v
                elif k == 'is_obsolete':
                    obsolete = True
                elif k == "is_a" and k in valid_rel:
                    s = v.split('!')[0].strip()
                    rel.append(s)
                elif k == "relationship" and v.startswith("part_of") and "part_of" in valid_rel:
                    s = v.split()[1].strip()
                    rel.append(s)

        # Last record
        if obsolete is False and namespace is not None:
            go_dict.setdefault(namespace, {})[term_id] = {'name': name,
                                                          'namespace': namespace,
                                                          'def': term_def,
                                                          'alt_id': alt_id,
                                                          'rel': rel}

    return go_dict


def gt_parser(gt_file, ontologies, ns_dict):
    """
    Parse ground truth file
    Discard terms not included in the ontology
    """

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

    gts = {}
    for ont in ontologies:
        if gt_dict.get(ont.namespace):
            terms = sorted([(v['index'], k, v['name']) for k, v in ont.go_terms.items()])
            matrix = np.zeros((len(gt_dict[ont.namespace]), ont.idxs), dtype='bool')

            ids = {}
            for i, p_id in enumerate(gt_dict[ont.namespace]):
                ids[p_id] = i
                for go_id in gt_dict[ont.namespace][p_id]:
                    matrix[i, ont.go_terms[go_id]['index']] = 1
            gts[ont.namespace] = GroundTruth(ids, matrix, terms=terms)
            logging.info('Ground truth: {}, proteins {}'.format(ont.namespace, i))

    return gts


def pred_parser(pred_file, ontologies, gts, ns_dict):
    """
    Parse a prediction file and returns a list of prediction objects, one for each namespace
    """
    pred_dict = {}
    with open(pred_file) as f:
        for line in f:
            line = line.strip().split()
            if line and len(line) > 2:
                p_id, go_id, prob = line[:3]
                if gts.get(ns_dict.get(go_id)) and p_id in gts[ns_dict[go_id]].ids:
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
            logging.info("Prediction: {}, {}, proteins {}".format(pred_file, ns, len(pred_dict[ns])))

    return predictions


def ia_parser(file):
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

from graph import Prediction, GroundTruth, propagate
import numpy as np
import logging
# import xml.etree.ElementTree as ET


def obo_parser(obo_file, valid_rel=("is_a", "part_of")):
    """
    Parse a OBO file and returns a list of ontologies, one for each namespace.
    Obsolete terms are excluded as well as external namespaces.
    """
    term_dict = {}
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
                        term_dict.setdefault(namespace, {})[term_id] = {'name': name,
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
            term_dict.setdefault(namespace, {})[term_id] = {'name': name,
                                                          'namespace': namespace,
                                                          'def': term_def,
                                                          'alt_id': alt_id,
                                                          'rel': rel}
    return term_dict


def gt_parser(gt_file, ontologies):
    """
    Parse ground truth file. Discard terms not included in the ontology.
    """
    gt_dict = {}
    with open(gt_file) as f:
        for line in f:
            line = line.strip().split()
            if line:
                p_id, term_id = line[:2]
                for ont in ontologies:
                    if term_id in ont.terms_dict:
                        gt_dict.setdefault(ont.namespace, {}).setdefault(p_id, []).append(term_id)
                        break

    gts = {}
    for ont in ontologies:
        if gt_dict.get(ont.namespace):
            matrix = np.zeros((len(gt_dict[ont.namespace]), ont.idxs), dtype='bool')
            ids = {}
            for i, p_id in enumerate(gt_dict[ont.namespace]):
                ids[p_id] = i
                for term_id in gt_dict[ont.namespace][p_id]:
                    matrix[i, ont.terms_dict[term_id]['index']] = 1
            logging.debug("gt matrix {} {} ".format(ont.namespace, matrix))
            propagate(matrix, ont, ont.order, mode='max')
            logging.debug("gt matrix propagated {} {} ".format(ont.namespace, matrix))
            gts[ont.namespace] = GroundTruth(ids, matrix, ont.namespace)
            logging.info('Ground truth: {}, proteins {}'.format(ont.namespace, len(ids)))

    return gts


def pred_parser(pred_file, ontologies, gts, prop_mode, max_terms=None):
    """
    Parse a prediction file and returns a list of prediction objects, one for each namespace.
    If a predicted is predicted multiple times for the same target, it stores the max.
    This is the slow step if the input file is huge, ca. 1 minute for 5GB input on SSD disk.

    """
    ids = {}
    matrix = {}
    ns_dict = {}  # {namespace: term}
    onts = {ont.namespace: ont for ont in ontologies}
    for ns in gts:
        matrix[ns] = np.zeros(gts[ns].matrix.shape, dtype='float')
        ids[ns] = {}
        for term in onts[ns].terms_dict:
            ns_dict[term] = ns

    with open(pred_file) as f:
        for line in f:
            line = line.strip().split()
            if line and len(line) > 2:
                p_id, term_id, prob = line[:3]
                ns = ns_dict.get(term_id)
                if ns in gts and p_id in gts[ns].ids:
                    i = gts[ns].ids[p_id]
                    if max_terms is None or np.count_nonzero(matrix[ns][i]) <= max_terms:
                        j = onts[ns].terms_dict.get(term_id)['index']
                        ids[ns][p_id] = i
                        matrix[ns][i, j] = max(matrix[ns][i, j], float(prob))

    predictions = []
    for ns in ids:
        if ids[ns]:
            logging.debug("pred matrix {} {} ".format(ns, matrix))
            propagate(matrix[ns], onts[ns], onts[ns].order, mode=prop_mode)
            logging.debug("pred matrix {} {} ".format(ns, matrix))

            predictions.append(Prediction(ids[ns], matrix[ns], len(ids[ns]), ns))
            logging.info("Prediction: {}, {}, proteins {}".format(pred_file, ns, len(ids[ns])))

    if not predictions:
        # raise Exception("Empty prediction, check format")
        logging.warning("Empty prediction! Check format or overlap with ground truth")

    return predictions


def ia_parser(file):
    ia_dict = {}
    with open(file) as f:
        for line in f:
            if line:
                term, ia = line.strip().split()
                ia_dict[term] = float(ia)
    return ia_dict


# def parse_sprot_xml(input_file, output_file):
#     """
#     Parse the Swiss-Prot XML annotation file
#     and write a TSV file with:
#     - accession
#     - GO term
#     - evidence code (ECO)
#     """
#     namespaces = {'uniprot': 'http://uniprot.org/uniprot'}
#     nsl = len(namespaces['uniprot']) + 2
#     with open(input_file) as f:
#         with open(output_file, 'w') as fout:
#             for event, elem in ET.iterparse(f, events=('start', 'end')):
#                 if event == 'end':
#                     if elem.tag[nsl:] == 'entry':
#                         acc = elem.find('uniprot:accession', namespaces).text
#                         for el in elem.iterfind('uniprot:dbReference', namespaces):
#                             if el.attrib['type'] == 'GO':
#                                 for at in el.iter():
#                                     if at.attrib['type'] == 'evidence':
#                                         fout.write('{}\t{}\t{}\n'.format(acc, el.attrib['id'], at.attrib['value']))
#                         elem.clear()
#     return

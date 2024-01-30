from cafaeval.graph import Graph, Prediction, GroundTruth, propagate
import numpy as np
import logging
logging.getLogger(__name__).addHandler(logging.NullHandler())
# import xml.etree.ElementTree as ET


def obo_parser(obo_file, valid_rel=("is_a", "part_of"), ia_file=None, orphans=True):
    """
    Parse a OBO file and returns a list of ontologies, one for each namespace.
    Obsolete terms are excluded as well as external namespaces.
    """
    # Parse the OBO file and creates a different graph for each namespace
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

    # Parse IA file
    ia_dict = None
    if ia_file is not None:
        ia_dict = ia_parser(ia_file)

    ontologies = {}
    for ns, ont_dict in term_dict.items():
        ontologies[ns] = Graph(ns, ont_dict, ia_dict, orphans)

    return ontologies


def gt_parser(gt_file, ontologies):
    """
    Parse ground truth file. Discard terms not included in the ontology.
    """
    gt_dict = {}
    replaced = {}
    with open(gt_file) as f:
        for line in f:
            line = line.strip().split()
            if line:
                p_id, term_id = line[:2]
                for ns in ontologies:
                    if term_id in ontologies[ns].terms_dict:
                        gt_dict.setdefault(ns, {}).setdefault(p_id, []).append(term_id)
                        break
                    # Replace alternative ids with canonical ids
                    elif term_id in ontologies[ns].terms_dict_alt:
                        for t_id in ontologies[ns].terms_dict_alt[term_id]:
                            gt_dict.setdefault(ns, {}).setdefault(p_id, []).append(t_id)
                            replaced.setdefault(ns, 0)
                            replaced[ns] += 1
                        break

    gts = {}
    for ns in ontologies:
        if gt_dict.get(ns):
            matrix = np.zeros((len(gt_dict[ns]), ontologies[ns].idxs), dtype='bool')
            ids = {}
            for i, p_id in enumerate(gt_dict[ns]):
                ids[p_id] = i
                for term_id in gt_dict[ns][p_id]:
                    matrix[i, ontologies[ns].terms_dict[term_id]['index']] = 1
            logging.debug("gt matrix {} {} ".format(ns, matrix))
            propagate(matrix, ontologies[ns], ontologies[ns].order, mode='max')
            logging.debug("gt matrix propagated {} {} ".format(ns, matrix))
            gts[ns] = GroundTruth(ids, matrix, ns)
            logging.info('Ground truth: {}, proteins {}, annotations {}, replaced alt. ids {}'.format(ns, len(ids),
                                                                                np.count_nonzero(matrix), replaced.get(ns, 0)))

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
    replaced = {}
    for ns in gts:
        matrix[ns] = np.zeros(gts[ns].matrix.shape, dtype='float')
        ids[ns] = {}
        for term in ontologies[ns].terms_dict:
            ns_dict[term] = ns
        for term in ontologies[ns].terms_dict_alt:
            ns_dict[term] = ns

    with (open(pred_file) as f):
        for line in f:
            line = line.strip().split()
            if line and len(line) > 2:
                p_id, term_id, prob = line[:3]
                ns = ns_dict.get(term_id)
                if ns in gts and p_id in gts[ns].ids:
                    # Get protein index
                    i = gts[ns].ids[p_id]
                    # Replace alternative ids with canonical ids
                    term_ids = [term_id]
                    if term_id in ontologies[ns].terms_dict_alt:
                        term_ids = ontologies[ns].terms_dict_alt[term_id]
                        replaced.setdefault(ns, 0)
                        replaced[ns] += len(term_ids)
                    for term_id in term_ids:
                        if max_terms is None or np.count_nonzero(matrix[ns][i]) <= max_terms:
                            j = ontologies[ns].terms_dict.get(term_id)['index']
                            ids[ns][p_id] = i
                            matrix[ns][i, j] = max(matrix[ns][i, j], float(prob))

    predictions = {}
    for ns in ids:
        if ids[ns]:
            logging.debug("pred matrix {} {} ".format(ns, matrix))
            propagate(matrix[ns], ontologies[ns], ontologies[ns].order, mode=prop_mode)
            logging.debug("pred matrix {} {} ".format(ns, matrix))

            predictions[ns] = Prediction(ids[ns], matrix[ns], ns)
            logging.info("Prediction: {}, {}, proteins {}, annotations {}, replaced alt. ids {}".format(pred_file, ns, len(ids[ns]),
                                                                                np.count_nonzero(matrix[ns]), replaced.get(ns, 0)))

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

import numpy as np
import copy
import logging


class Graph:
    """
    Ontology class. One ontology == one namespace
    DAG is the adjacence matrix (sparse) which represent a Directed Acyclic Graph where
    DAG(i,j) == 1 means that the go term i is_a (or is part_of) j
    Parents that are in a different namespace are discarded
    """
    def __init__(self, namespace, go_dict, ia_dict=None, orphans=False):
        """
        go_dict = {term: {name: , namespace: , def: , alt_id: , rel:}}
        """
        self.namespace = namespace
        self.dag = []  # [[], ...] terms (rows, axis 0) x parents (columns, axis 1)
        self.go_terms = {}  # {term: {index: , name: , namespace: , def: }  used to assign term indexes in the gt
        self.go_list = []  # [{id: term, name:, namespace: , def:, adg: [], children: []}, ...]
        self.idxs = None  # Number of terms
        self.order = None
        self.toi = None
        self.ia = None

        rel_list = []
        for self.idxs, (go_id, term) in enumerate(go_dict.items()):
            rel_list.extend([[go_id, rel, term['namespace']] for rel in term['rel']])
            self.go_list.append({'id': go_id, 'name': term['name'], 'namespace': namespace, 'def': term['def'],
                                 'adj': [], 'children': []})
            self.go_terms[go_id] = {'index': self.idxs, 'name': term['name'], 'namespace': namespace, 'def': term['def']}
            for a_id in term['alt_id']:
                self.go_terms[a_id] = copy.copy(self.go_terms[go_id])
        self.idxs += 1

        self.dag = np.zeros((self.idxs, self.idxs), dtype='bool')

        # id1 term (row, axis 0), id2 direct parent (column, axis 1)
        for id1, id2, ns in rel_list:
            if self.go_terms.get(id2):
                i = self.go_terms[id1]['index']
                j = self.go_terms[id2]['index']
                self.dag[i, j] = 1
                self.go_list[i]['adj'].append(j)
                self.go_list[j]['children'].append(i)
            else:
                logging.debug('Skipping branch to external namespace: {}'.format(id2))

        self.order = top_sort(self)

        if orphans:
            self.toi = np.fill(self.dag.shape[1], 1)  # All terms, also those without parents
        else:
            self.toi = np.where(self.dag.sum(axis=1) > 0)[0]  # Only terms with parents

        if ia_dict is not None:
            self.set_ia(ia_dict)

        return

    def set_ia(self, ia_dict):
        self.ia = np.zeros(self.idxs, dtype='float')
        for go_id in self.go_terms:
            if ia_dict.get(go_id):
                self.ia[self.go_terms[go_id]['index']] = ia_dict.get(go_id)
            else:
                logging.debug('Missing IA for term: {}'.format(go_id))
        # Convert inf to zero
        np.nan_to_num(self.ia, copy=False, nan=0, posinf=0, neginf=0)


class Prediction:
    """
    The score matrix contains the scores given by the predictor for every node of the ontology
    """
    def __init__(self, ids, matrix, idx, namespace=None):
        self.ids = ids
        self.matrix = matrix  # scores
        self.next_idx = idx
        self.n_pred_seq = idx + 1
        self.namespace = namespace

    def __str__(self):
        return "\n".join(["{}\t{}\t{}".format(index, self.matrix[index], self.namespace) for index, _id in enumerate(self.ids)])


class GroundTruth:
    def __init__(self, ids, matrix, terms=None):
        self.ids = ids
        self.matrix = matrix
        self.terms = terms


def top_sort(ont):
    """
    Takes a sparse matrix representing a DAG and returns an array with the node indexes in topological order
    https://en.wikipedia.org/wiki/Topological_sorting
    """
    indexes = []
    visited = 0
    (rows, cols) = ont.dag.shape
    
    # create a vector containing the in-degree of each node 
    in_degree = []
    for c in range(0, cols):
        in_degree.append(ont.dag[:, c].sum())
    # find the nodes with in-degree 0 and add them to the queue
    queue = [index for index, value in enumerate(in_degree) if value == 0]

    # for each element of the queue increment visits, add them to the list of ordered nodes and decrease the in-degree
    # of the neighbor nodes and add them to the queue if they reach in-degree == 0
    while queue:
        visited += 1
        idx = queue.pop(0)
        indexes.append(idx)
        in_degree[idx] -= 1
        l = ont.go_list[idx]['adj']
        if len(l) > 0:
            for j in l:
                in_degree[j] -= 1
                if in_degree[j] == 0:
                    queue.append(j)
    # if visited is equal to the number of nodes in the graph then the sorting is complete 
    # otherwise the graph can't be sorted with topological order
    if visited == rows:
        return indexes
    else:
        raise Exception("The sparse matrix doesn't represent an acyclic graph")


def propagate(gt, ont, order):
    """
    Takes the score matrix (proteins x terms) and
    a sparse matrix representing the DAG and returns a new score matrix
    with terms propagated upwards
    """
    if gt.matrix.shape[0] == 0:
        raise Exception("Empty ground truth")

    deepest = np.where(np.sum(gt.matrix[:, order], axis=0) > 0)[0][0]
    if deepest.size == 0:
        raise Exception("The matrix is empty")
    order = np.delete(order, [range(0, deepest)])

    for i in order:
        current = np.where(ont.dag[:, i] != 0)[0]
        if current.size > 0:
            cols = np.concatenate((current, [i]))
            gt.matrix[:, i] = gt.matrix[:, cols].max(axis=1)

    return order





import numpy as np
import copy
import logging
from evaluation import get_toi_idx


class Graph:
    # TODO term GO:0000030 has parent in both BPO and MFO, think about how to fix it
    """
    Ontology class. One ontology == one namespace
    DAG is the adjacence matrix (sparse) which represent a Directed Acyclic Graph where
    DAG(i,j) = 1 means that the go term i is_a (or is part_of) j
    """
    def __init__(self, namespace, go_dict, ia_dict=None):
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
        self.ia_array = None

        idxs = 0
        rel_list = []
        for go_id, term in go_dict.items():
            rel_list.extend([[go_id, rel, term['namespace']] for rel in term['rel']])
            self.go_list.append({'id': go_id, 'name': term['name'], 'namespace': namespace, 'def': term['def'],
                                 'adj': [], 'children': []})
            self.go_terms[go_id] = {'index': idxs, 'name': term['name'], 'namespace': namespace, 'def': term['def']}
            for a_id in term['alt_id']:
                self.go_terms[a_id] = copy.copy(self.go_terms[go_id])
            idxs += 1

        self.idxs = idxs
        self.dag = np.zeros((idxs, idxs), dtype='bool')

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
        self.toi = get_toi_idx(self.dag)

        if ia_dict is not None:
            self.set_ia_array(ia_dict)

        return

    def set_ia_array(self, ic_dict):
        n = len(self.go_list)
        self.ia_array = np.zeros(n, dtype='float')
        for i in range(0, n):
            self.ia_array[i] = ic_dict.get(self.go_list[i]['id'], 0.0)


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


# takes a sparse matrix representing a DAG and returns an array containing the node indexes in topological order
def top_sort(ont):
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


# takes the score matrix (proteins on rows and terms on columns) and
# a sparse matrix representing the DAG and returns a new score matrix
# with terms propagated upwards
def propagate(gt, ont, order=None):
    if gt.matrix.shape[0] == 0:
        raise Exception("Empty ground truth")

    if order is None:
        order = top_sort(ont)
    
    deepest = np.where(np.sum(gt.matrix[:, order], axis=0) > 0)[0][0]
    if deepest.size == 0:
        raise Exception("The matrix is empty")
    order = np.delete(order, [range(0, deepest)])

    for i in order:
        C = np.where(ont.dag[:, i] != 0)[0]
        if C.size > 0:
            cols = np.concatenate((C, [i]))
            gt.matrix[:, i] = gt.matrix[:, cols].max(axis=1)
    return order





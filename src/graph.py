import numpy as np


class Graph:
    """
    Class used to contain an ontology
    DAG is the adjacence matrix (sparse) which represent a Directed Acyclic Graph where
    DAG(i,j) = 1 means that the go term i is_a (or is part_of) j
    """
    def __init__(self, dag, go_terms, ontology_type, go_list, ia_dict=None):
        self.dag = dag  # [[], ...]  terms X terms matrix
        self.go_terms = go_terms  # {term: {index: , name: , namespace: , def: }
        self.ontology_type = ontology_type  # namespace
        self.go_list = go_list  # [{id: term, name:, namespace: , def:, adg: [], children: []}, ...]
        self.ia_array = None
        if ia_dict is not None:
            self.set_ia_array(ia_dict)

    def get_id(self, index):
        return self.go_terms[index]['id']

    def get_namespace(self):
        return self.ontology_type

    def set_ia_array(self, ic_dict):
        n = len(self.go_list)
        self.ia_array = np.zeros(n, dtype='float')
        for i in range(0, n):
            self.ia_array[i] = ic_dict.get(self.go_list[i]['id'], 0.0)


class Prediction:
    """
    The score matrix contains the scores given by the predictor for every node of the ontology
    """
    def __init__(self, ids, scores, idx, namespace=None):
        self.ids = ids
        self.matrix = scores
        self.next_idx = idx
        self.n_pred_seq = idx + 1
        self.namespace = namespace

    def __str__(self):
        return "\n".join(["{}\t{}\t{}".format(index, self.matrix[index], self.namespace) for index, _id in enumerate(self.ids)])


class GroundTruth:
    def __init__(self, p_ids, gt_matrix, terms=None):
        self.ids = p_ids
        self.matrix = gt_matrix
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
        return
    if order is None:
        order = top_sort(ont)
    
    deepest = np.where(np.sum(gt.matrix[:, order], axis=0) > 0)[0][0]
    if deepest.size == 0:
        raise Exception("The matrix is empty")
    order = np.delete(order, [range(0,deepest)])

    for i in order:
        C = np.where(ont.dag[:, i] != 0)[0]
        if C.size > 0:
            cols = np.concatenate((C, [i]))
            gt.matrix[:, i] = gt.matrix[:, cols].max(axis=1)
    return order


    



"""scores = np.array([[0.5,0,0], [0,0.5,0], [0,0,0.5]])
dag = np.array([[0,0,0], [1,0,0], [0,1,0]])

print("scores \n", scores)
prop, order = propagate_pred(scores,dag)
print(order)
print("propagated scores \n", prop)
print("propagated scores \n", prop[:,order])"""

"""ont = np.array([[0,1,1], [0,0,0], [0,0,0]]) #problema girando le matrici
print("order", top_sort(ont))
print("propagated scores \n", propagate_pred(scores,ont))"""






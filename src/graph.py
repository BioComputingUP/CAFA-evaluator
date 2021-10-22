import numpy as np
import time

# DAG is the adjacence matrix (sparse) which represent a Directed Acyclic Graph where DAG(i,j) = 1
# means that the go term i is part of j

## The score matrix contains the scores given by the predictor for every node of the ontology

# class used to contain an ontology 
class Graph:
    def __init__(self, dag, go_terms, ontology_type, go_list):
        self.dag = dag 
        self.go_terms = go_terms
        self.ontology_type = ontology_type
        self.go_list = go_list

    def get_id(index):
        return go_terms[index]['id']

    def get_namespace(self):
        return self.ontology_type

class Prediction:
    def __init__(self, ids, scores, idx, namespace=None):
        self.ids = ids
        self.scores = scores
        self.next_idx = idx
        self.n_pred_seq = idx + 1
        self.namespace = namespace
        #self.unused_scores = unused_scores


class Ground_truth:
    def __init__(self, p_ids, gt_matrix):
        self.ids = p_ids
        self.matrix = gt_matrix



# takes a sparse matrix representing a DAG and returns an array containing the node indexes in topological order
def top_sort(ont):
    indexes = []
    visited = 0
    (rows, cols) = ont.dag.shape
    
    # create a vector containing the in-degree of each node 
    in_degree = []
    for c in range(0,cols):
        in_degree.append(ont.dag[:,c].sum())
    #find the nodes with in-degree 0 and add them to the queue 
    queue = [index for index, value in enumerate(in_degree) if value == 0]

    # for each element of the queue increment visits, add them to the list of ordered nodes and decrease the in-degree
    #of the neighbor nodes and add them to the queue if they reach in-degree == 0
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


# takes the score matrix (protein on rows and go terms on columns) and a sparse matrix representing a DAG
# and returns a the score matrix with ptedictions propagated upwards
def propagate_pred(pred, ont, order=None):
    if pred.scores.shape[0] == 0:
        return
    # get the go term indexes in topological order
    if order is None:
        order = top_sort(ont)
    # find the the deepest node with a score assigned to it 
    #print(pred.scores.sum(axis=1).sum())
    deep_begin = time.time()
    deepest = np.where(np.sum(pred.scores[:,order], axis=0) > 0)[0]
    if deepest.size != 0:
        #raise Exception("The score matrix contains no predictions") 
        order = np.delete(order, [range(0,deepest[0])]) 
    deep_end = time.time()
    #else:
        #print(np.where(np.sum(pred.scores[:,order], axis=0) > 0))
        #print(pred.ids)

    #deepest = np.where(np.sum(scores, axis=0) > 0)[0]
    #print(deepest)
    # for each element in order take find the childes and update the scores on its column with the max score between
    # its score and the childs score 
    prop_begin = time.time()
    for i in order: #range(0,dag.shape[1]): #range(0,len(order))
        # children search
        #C = np.where(ont.dag[:,i] != 0)[0]
        C = ont.go_list[i]['children']
        if len(C) > 0:
            cols = np.concatenate((C, [i]))
            pred.scores[:,i] = pred.scores[:,cols].max(axis=1)
    prop_end = time.time()
    return order

def propagate_gt(gt, ont, order=None):
    if gt.matrix.shape[0] == 0:
        return
    if order is None:
        order = top_sort(ont)
    
    deepest = np.where(np.sum(gt.matrix[:,order], axis=0) > 0)[0][0]
    if deepest.size == 0:
        raise Exception("The ground truth matrix is empty")
    order = np.delete(order, [range(0,deepest)])

    for i in order:
        C = np.where(ont.dag[:,i] != 0)[0]
        if C.size > 0:
            cols = np.concatenate((C, [i]))
            gt.matrix[:,i] = gt.matrix[:,cols].max(axis=1)
    return order

def ic_to_array(ic_dict, ont):
    n = len(ont.go_list)
    ic_arr = np.zeros(n, dtype='float')
    for i in range(0,n):
        ic = ic_dict.get(ont.go_list[i]['id'])
        if ic is not None:
            ic_arr[i] = ic
    return ic_arr
    



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






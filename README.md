# CAFA-evaluator

Calculate precision-recall curves and F-max as in the 
[Critical Assessment of protein Function Annotation (CAFA)](https://www.biofunctionprediction.org/cafa/)

CAFA-evaluator is generic and works with any ontology. For example it is used to evaluate 
intrinsic disorder function predictions against [DisProt](https://disprot.org/) annotations.

## Required packages

- numpy
- logging
- matplotlib

## HowTo

To run CAFA-evaluator simply call **main.py** in the src/ folder with the following arguments:

Mandatory
* *Ontology file* containing the structure of the terms in OBO format

* *Prediction folder* with prediction files. Sub-folders are processed recursively and the sub-folder name is used as prefix

* *Ground truth folder* with ground truth files, one for each namespace (sub-ontology). 
The ontology namespace has to be hardcoded in the filename, e.g. gt_*bpo*_.txt or *disorderfunction*.tsv

Optional

* *Output folder*
* *Benchmark folder* with benchmark files, one for each namespace (sub-ontology). 
The evaluation is performed considering only the targets listed in these files. 
The ontology namespace has to be hardcoded in the filename, e.g. bm_*bpo*_.txt or *disorderfunction*.tsv



### Namespaces (sub-ontologies)

At the moment in order to recognize the namespaces you should 
modify this dictionary which is hardcoded in *main.py*.

The key is the string matched in the ground truth and benchmark files, 
while values correspond to the namespace field in the OBO file. 

~~~Python3
namespaces = {"bpo": "biological_process", "cco": "cellular_component", "mfo": "molecular_function",
 "disorderfunction": "Disorder function", "interactionpartner": "Interaction partner",
 "structuralstate": "Structural state", "structuraltransition": "Structural transition"}
~~~



### Output

* Precision-recall graph for each ontology (<namespace>.png)
* Text file with precision, recall and Fmax values for each namespace (results.tsv)
* Log file (info.log)

### Prediction format

Target_ID   Term_ID Score

~~~txt
T98230000001    GO:0000159      0.39
T98230000001    GO:0032991      0.39
T98230000001    GO:1902494      0.39
...
~~~

### Ground truth format

Target_ID   Term_ID	

~~~
T100900005305   GO:0033234
T100900010085   GO:0006468
T100900010085   GO:0046777
...
~~~


### Benchmark format

Path to the folder containing the files, one for each ontology, in which are listed all the proteins that will be actually used for evaluation (the '/' after the name of the folder must be ommitted). like in the ground truth case each file has to have the namespace in the name of the file i.e bpo_benchmark.txt. A filter file should just be a list of protein Ids like this:

~~~
T446890002807
T446890001819
T446890003376
...
~~~



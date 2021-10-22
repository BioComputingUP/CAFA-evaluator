# CAFA-evaluator

## Required packages

- Numpy
- Logging
- matplotlib

## HowTo

To run the evaluation on some existing predictions run **main.py** in src with the following arguments:

### Ontology file

File .obo that contains the structure of the ontology.

### Prediction Folder

The path to the folder containing the predictions (the '/' after the name of the folder must be ommitted). The predictions have to be contained in a folder nested in the given path to allow multiple predictior comparison, for example if we give the path to the folder *Prediction* as parameter the structure has to be:

~~~
Predictions    
│
└───PredictorA
│   │   predictions1.txt
│   │   predictions2.txt
│   │	...
│   
└───PredictorB
    │   predictions1.txt
    │   predictions2.txt
    |   ...
~~~

The predictions can be contained in a single file but they can be split up in more files in case of memory problems. The only requirement is that all the predictions for a specific protein must be contained in a single file. The format of prediction files has to be:

Cafa_id	Predicted_Go_Term	Score

~~~txt
T98230000001    GO:0000159      0.39
T98230000001    GO:0032991      0.39
T98230000001    GO:1902494      0.39
...
~~~



### Ground Truth Folder

Path to the folder containing the ground truth files one for each ontology used in the evaluation (the '/' after the name of the folder must be ommitted). In order for the program to pick up this files their name has to contain the namespace, for example the ground truth file for the biological process should be named *bpo_gt.txt*. The format of a ground truth file should be:

Cafa_Id	Go_Term	

~~~
T100900005305   GO:0033234
T100900010085   GO:0006468
T100900010085   GO:0046777
...
~~~

### (optional) output path

Path to the folder where the results will be written.

### (optional) gt_filters

Path to the folder containing the files, one for each ontology, in which are listed all the proteins that will be actually used for evaluation (the '/' after the name of the folder must be ommitted). like in the ground truth case each file has to have the namespace in the name of the file i.e bpo_benchmark.txt. A filter file should just be a list of protein Ids like this:

~~~
T446890002807
T446890001819
T446890003376
...
~~~

## Allowed namespaces

At the moment the namespaces code are hard-coded in main like this:

~~~Python3
namespaces = {"bpo": "biological_process", "cco": "cellular_component", "mfo": "molecular_function",
 "disorderfunction": "Disorder function", "interactionpartner": "Interaction partner",
 "structuralstate": "Structural state", "structuraltransition": "Structural transition"}
~~~

Where the key is the code that will needs to be in the ground-truth and filter file meanwhile the values correspond to the namespace field in the .obo file. If more ontologies need to be used just modify this dictionary. 	

## Outputs

- A precision-recall graph for each ontology 
- A results.tsv file containing the precision, recall and fmax for each ontology
- A log file containing some information about terms that are not used 


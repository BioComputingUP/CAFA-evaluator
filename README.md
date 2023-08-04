# CAFA-evaluator

CAFA-evaluator is a Python program designed to evaluate the performance of prediction methods on targets 
with hierarchical concept dependencies.
It generalizes multi-label evaluation to modern ontologies where the prediction targets are drawn 
from a directed acyclic graph and achieves high efficiency by leveraging matrix computation and topological sorting.

The code replicates the 
[Critical Assessment of protein Function Annotation (CAFA)](https://www.biofunctionprediction.org/cafa/) benchmarking, 
which evaluates predictions of the consistent subgraphs in [Gene Ontology](http://geneontology.org/). 
The package also contains a **Jupyter Notebook** to generate precision-recall and remaining uncertainty–misinformation curves.
CAFA-evaluator implementation was inspired by the original Matlab code used in CAFA2 assessment and available at 
[https://github.com/yuxjiang/CAFA2](https://github.com/yuxjiang/CAFA2).

## Installation

From GitHub:

    git clone https://github.com/BioComputingUP/CAFA-evaluator.git  
    pip install .

From PyPI:

    pip install cafaeval

## Usage

Example input files are provided inside the `data/example` folder. If you have installed the software using pip,
the folder is located in the Python site-packages directory.

**Note**: Weighted scores are generated only if the *Information Accretion* file is provided.


### Library

The `cafa_eval` function is the main entry point of the package. It accepts the following arguments:

* Ontology file in OBO format

* Prediction folder contain prediction files. Files can be organized into sub-folders, 
sub-folders are processed recursively and the sub-folder name is used as prefix for the method name

* Ground truth file containing targets and associated ontology terms

Optional arguments are:

* Information accretion file
* Flag to exclude orphan terms, e.g. the root(s)
* Normalization strategy
* Propagation strategy 
* Max number of terms to consider for each protein and namespace
* Threshold step size
* Max number of threads to use

```pycon
>>> import cafaeval
>>> from cafaeval.evaluation import cafa_eval
>>> cafa_eval("IDPO_disorder_function.obo", "predictions", "ground_truth.tsv")
(                                        cov        pr        rc         f
filename   ns                tau                                         
pred_5.tsv disorder_function 0.01  1.000000  0.292532  0.959623  0.448380
                             0.02  1.000000  0.292532  0.959623  0.448380
                             0.03  1.000000  0.292695  0.959623  0.448571
                             0.04  1.000000  0.292695  0.959623  0.448571
                             0.05  1.000000  0.292695  0.959623  0.448571
...                                     ...       ...       ...       ...
pred_1.tsv disorder_function 0.41  0.005952  0.250000  0.001984  0.003937
                             0.42  0.005952  0.250000  0.001984  0.003937
                             0.43  0.005952  0.250000  0.001984  0.003937
                             0.44  0.005952  0.250000  0.001984  0.003937
                             0.45  0.005952  0.333333  0.001984  0.003945

[352 rows x 4 columns], {'f':                                         cov        pr        rc         f  max_cov
filename   ns                tau                                                  
pred_1.tsv disorder_function 0.04  0.988095  0.466566  0.579960  0.517120      1.0
pred_2.tsv disorder_function 0.84  0.970238  0.504499  0.579960  0.539604      1.0
pred_3.tsv disorder_function 0.89  1.000000  0.638889  0.701290  0.668637      1.0
pred_4.tsv disorder_function 0.06  1.000000  0.777778  0.774504  0.776137      1.0
pred_5.tsv disorder_function 0.38  0.994048  0.596671  0.776389  0.674768      1.0})
```
Results can be saved to file using the `save_results` function:

```pycon
>>> import cafaeval
>>> from cafaeval.evaluation import cafa_eval, write_results
>>> res = cafa_eval("IDPO_disorder_function.obo", "predictions", "ground_truth.tsv")
>>> write_results(*res)

```

The output of `cafa_eval` is a tuple containing:
* A pandas DataFrame with the evaluation results, one row per prediction file, namespace and threshold.
* A dictionary with the best scores (max F-measure, max Weighted F-measure, min Semantic similarity). For each 
score the dictionary contain a pandas DataFrame with one row per prediction file, namespace and threshold. 

The `write_results` function generated the following files:
* *evaluation_all.tsv* corresponding to the first object returned by the `cafa_eval` function.
* *evaluation_best_< metric >.tsv* corresponding to the second object returned by the `cafa_eval` function. 
A different file for each metric is created.

### Command line

When executed from the command line the script logs information about the calculation in the console (standard error) and
will create a folder named `results` containing the evaluation results. 
A different folder can be specified using the `-out_dir` option. 

```bashcon
cafaeval IDPO_disorder_function.obo predictions ground_truth.tsv 
```


## Input files
**Prediction file** - Tab separated file with the target ID, term ID and score columns.

~~~txt
T_1	IDPO:00501	0.06
T_1	IDPO:00506	0.05
T_1	IDPO:00507	0.03
T_2	IDPO:00501	0.04
T_2	IDPO:00506	0.02
...
~~~

**Ground truth file** - Tab separated file with the target ID and term ID. 
Additional columns are discarded.
~~~
T_1	IDPO:00024
T_2	IDPO:00506
T_3	IDPO:00502
T_4	IDPO:00025
...
~~~

**Information accretion (optional)** - If not provided the weighted and S statistics are not generated.
Information accretion (IA) can be calculated as described in
[Wyatt and Radivojac, Bioinformatics, 2013](https://pubmed.ncbi.nlm.nih.gov/23813009/).

```
IDPO:00024  6.32
IDPO:00506  12.04
IDPO:00502  1.34
IDPO:00025  0.56
...
```


## Plotting

Plots are generated by running all the cells in the `plot.ipynb` Jupyter Notebook.
In order to generate the figures you need to manually modify the first cell of the notebook which 
contains a few parameters: 
* the path to the input file generated in the evaluation`evaluation_all.tsv`, see [Usage](#usage) above.
* the output folder
* the metric (F-score, weighted F-score, S measure).
* an optional file including information about methods aliases and groups. If provided, the results
are presented selecting only one method per group. The file should look like:
```
filename	group	label
pred_1.tsv	BioLab	BioLab_model_1
pred_2.tsv	BioLab	BioLab_model_2
pred_3.tsv	JohnLab	JohnLab_model_1
```

The notebook generates the following files:

* *fig_< metric >_< name_space >.png* A figure for each namespace in the dataframe and the selected metric. 
The notebook generates one metric at the time, you have to modify the input cell to generate the plots for a different metric

* fig_< metric >.tsv* - A file with the data points for the metric curves. One curve for each method.
```
group	label	ns	tau	cov	wrc	wpr	wf
INGA	INGA_1	biological_process	0.010	0.993	0.557	0.094	0.160
INGA	INGA_1	biological_process	0.020	0.993	0.555	0.094	0.161
INGA	INGA_1	biological_process	0.030	0.993	0.552	0.095	0.162
INGA	INGA_1	biological_process	0.040	0.993	0.551	0.095	0.163
```

## Algorithm information 

#### Notes
* The word `aspect`, `namespace` and `sub-ontology` are used interchangeably in the following documentation.
* In order to replicate CAFA results, you can simply adapt the input files. 
  - *No/partial knowledge* can be reproduced by filtering/splitting the **ground truth file** 
  - In order to exclude specific terms from the analyses, 
e.g. generic "binding" terms, you can directly modify the input **ontology file** 

#### Input filtering

Prediction files are filtered considering only those targets included in the ground truth and 
only those terms included in the ontology file. 
If the ground truth contains only annotations from one aspect (e.g. "molecular function"), 
the evaluation is provided only for that aspect.
The `-max_terms` parameter determines the maximum number of terms that will be considered for each target and namespace. 
Parsing stops when the target limit for every namespace is reached. The score is not checked, 
meaning that terms are not sorted before the check, and the check is performed before propagation.

The ontology is processed with an internal parser that accepts only the OBO format. 
The following rules are applied: 
  - Obsolete terms are always excluded
  - Only "is_a" and "part_of" relationships are considered
  - Cross-aspect (cross-namespace) relationships are always discarded
  - Alternative term IDs are automatically mapped to the main ID

When information accretion is provided, terms which are not available in the accretion file are 
removed from the ontology.

#### Terms propagation

Both the predictions and the ground truth annotations are always propagated up to the ontology root(s). 
Two strategies are available: i) prediction scores are propagated without overwriting the scores 
assigned to the parents; ii) scores are propagated considering always the max.

#### Memory usage

- The algorithm stores in memory a Numpy boolean *N x M* array 
(N = number of ground truth targets; M = ontology terms of a single aspect)
for each aspect in the ground truth file.

- An array of the same size (rows &le; N), but containing floats (the prediction scores) instead of booleans, 
is stored for each prediction file. Prediction files are processed one by one and the matrix gets reassigned.




## CAFA5 challenge (Kaggle)
Owing to its reliability and accuracy, the organizers have selected CAFA-evaluator as the official evaluation software
in the [CAFA5 Kaggle](https://www.kaggle.com/competitions/cafa-5-protein-function-prediction) competition. 
In Kaggle the software is executed with the following command:

    cafaeval go-basic.obo prediction_dir test_terms.tsv -ia IA.txt -prop fill -norm cafa -th_step 0.001 -max_terms 500

In the example above the method prediction file should be inside the `prediction_dir` folder and 
evaluated against the `test_terms.tsv` file (not available to participants) containing the ground truth.



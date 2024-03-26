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

Visit the [CAFA-evaluator Wiki](https://github.com/BioComputingUP/CAFA-evaluator/wiki) for more information about the algorithm.

## Citation
[CAFA-evaluator: A Python Tool for Benchmarking Ontological Classification Methods](https://doi.org/10.1093/bioadv/vbae043)  
*D Piovesan, D Zago, P Joshi, MC De Paolis Kaluza, M Mehdiabadi, R Ramola, AM Monzon, W Reade, I Friedberg, P Radivojac, SCE Tosatto*  
**Bioinformatics Advances (2024)** - DOI: [10.1093/bioadv/vbae043](https://doi.org/10.1093/bioadv/vbae043)

## Installation

You can use the tool simply cloning this repo and running `__main__.py` Python script. 
which implements a command line interface.
To avoid import errors you need to add the source root folder 
to the `PYTHONPATH` environment variable. For example:

    export PYTHONPATH=/path/to/CAFA-evaluator/src:$PYTHONPATH

Alternatively, you can install the package with Pip. No need to export any variable
in this case.

From GitHub:

    git clone https://github.com/BioComputingUP/CAFA-evaluator.git  
    pip install .

From PyPI:

    pip install cafaeval

## Usage

The program can be executing the command line interface or as a library.
Both the command line and the library accept the following positional arguments:

* **Ontology file** in OBO format

* **Prediction folder** contain prediction files. Files can be organized into sub-folders, 
sub-folders are processed recursively and the sub-folder name is used as prefix for the method name

* **Ground truth file** containing targets and associated ontology terms

Example input files are provided inside the `data/example` folder. If you have installed the software using pip,
the folder is located in the Python `site-packages` directory.

### Command line

When executed from the command line the script logs information about the calculation in the console (standard error) and
will create a folder named `results` containing the evaluation results. 
A different folder can be specified using the `-out_dir` option. 

The following options are available:
When installed with pip the script is available as `cafaeval` command.

```bashcon
cafaeval ontology_file prediction_folder ground_truth_file 
```
If you simply cloned the repository:
```bashcon
python3 /path/to/CAFA-evaluator/src/cafaeval/__main__.py ontology_file prediction_folder ground_truth_file
```

### Library

The `cafa_eval` function is the main entry point of the package. 
Below is reported an example using the example files provided in the `data/example` folder.


```pycon
>>> import cafaeval
>>> from cafaeval.evaluation import cafa_eval
>>> cafa_eval("IDPO_disorder_function.obo", "predictions", "ground_truth.tsv")
(                                       n        tp        fp        fn        pr  ...         f         s  pr_micro  rc_micro   f_micro
filename   ns                tau                                                  ...                                                  
pred_5.tsv disorder_function 0.01  168.0  2.928571  7.083333  0.130952  0.292532  ...  0.448380  7.084544  0.292509  0.957198  0.448087
                             0.02  168.0  2.928571  7.083333  0.130952  0.292532  ...  0.448380  7.084544  0.292509  0.957198  0.448087
                             0.03  168.0  2.928571  7.077381  0.130952  0.292695  ...  0.448571  7.078592  0.292683  0.957198  0.448292
                             0.04  168.0  2.928571  7.077381  0.130952  0.292695  ...  0.448571  7.078592  0.292683  0.957198  0.448292
                             0.05  168.0  2.928571  7.077381  0.130952  0.292695  ...  0.448571  7.078592  0.292683  0.957198  0.448292
...                                  ...       ...       ...       ...       ...  ...       ...       ...       ...       ...       ...
pred_1.tsv disorder_function 0.41    1.0  0.005952  0.017857  3.053571  0.250000  ...  0.003937  3.053624  0.250000  0.001946  0.003861
                             0.42    1.0  0.005952  0.017857  3.053571  0.250000  ...  0.003937  3.053624  0.250000  0.001946  0.003861
                             0.43    1.0  0.005952  0.017857  3.053571  0.250000  ...  0.003937  3.053624  0.250000  0.001946  0.003861
                             0.44    1.0  0.005952  0.017857  3.053571  0.250000  ...  0.003937  3.053624  0.250000  0.001946  0.003861
                             0.45    1.0  0.005952  0.011905  3.053571  0.333333  ...  0.003945  3.053595  0.333333  0.001946  0.003868

[352 rows x 14 columns], {'f':                                        n        tp        fp        fn        pr  ...         s  pr_micro  rc_micro   f_micro  cov_max
filename   ns                tau                                                  ...                                                 
pred_1.tsv disorder_function 0.04  166.0  1.744048  2.095238  1.315476  0.466566  ...  2.473964  0.454264  0.570039  0.505608      1.0
pred_2.tsv disorder_function 0.84  163.0  1.755952  1.845238  1.303571  0.504499  ...  2.259248  0.487603  0.573930  0.527256      1.0
pred_3.tsv disorder_function 0.89  168.0  2.113095  1.601190  0.946429  0.638889  ...  1.859983  0.568910  0.690661  0.623902      1.0
pred_4.tsv disorder_function 0.06  168.0  2.333333  0.666667  0.726190  0.777778  ...  0.985798  0.777778  0.762646  0.770138      1.0
pred_5.tsv disorder_function 0.38  167.0  2.345238  1.952381  0.714286  0.596671  ...  2.078941  0.545706  0.766537  0.637540      1.0

[5 rows x 15 columns], 's':                                        n        tp        fp        fn        pr  ...         s  pr_micro  rc_micro   f_micro  cov_max
filename   ns                tau                                                  ...                                                 
pred_1.tsv disorder_function 0.06  102.0  1.047619  0.297619  2.011905  0.816667  ...  2.033799  0.778761  0.342412  0.475676      1.0
pred_2.tsv disorder_function 0.91  124.0  1.196429  0.839286  1.863095  0.638978  ...  2.043410  0.587719  0.391051  0.469626      1.0
pred_3.tsv disorder_function 0.89  168.0  2.113095  1.601190  0.946429  0.638889  ...  1.859983  0.568910  0.690661  0.623902      1.0
pred_4.tsv disorder_function 0.06  168.0  2.333333  0.666667  0.726190  0.777778  ...  0.985798  0.777778  0.762646  0.770138      1.0
pred_5.tsv disorder_function 0.42  156.0  1.642857  0.714286  1.416667  0.776221  ...  1.586552  0.696970  0.536965  0.606593      1.0

[5 rows x 15 columns], 'f_micro':                                        n        tp        fp        fn        pr  ...         s  pr_micro  rc_micro   f_micro  cov_max
filename   ns                tau                                                  ...                                                 
pred_1.tsv disorder_function 0.04  166.0  1.744048  2.095238  1.315476  0.466566  ...  2.473964  0.454264  0.570039  0.505608      1.0
pred_2.tsv disorder_function 0.84  163.0  1.755952  1.845238  1.303571  0.504499  ...  2.259248  0.487603  0.573930  0.527256      1.0
pred_3.tsv disorder_function 0.89  168.0  2.113095  1.601190  0.946429  0.638889  ...  1.859983  0.568910  0.690661  0.623902      1.0
pred_4.tsv disorder_function 0.06  168.0  2.333333  0.666667  0.726190  0.777778  ...  0.985798  0.777778  0.762646  0.770138      1.0
pred_5.tsv disorder_function 0.38  167.0  2.345238  1.952381  0.714286  0.596671  ...  2.078941  0.545706  0.766537  0.637540      1.0

[5 rows x 15 columns]})
```

The output of `cafa_eval` is a tuple containing:
* A pandas DataFrame with the evaluation results, one row per prediction file, namespace and threshold.
* A dictionary with the best scores (max F-measure, max Weighted F-measure, min Semantic similarity). For each 
score the dictionary contain a pandas DataFrame with one row per prediction file, namespace and threshold. 

The `write_results` function generates the output files.

```pycon
>>> import cafaeval
>>> from cafaeval.evaluation import cafa_eval, write_results
>>> res = cafa_eval("IDPO_disorder_function.obo", "predictions", "ground_truth.tsv")
>>> write_results(*res)
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

**Information accretion file (optional)** - If not provided, the weighted and S statistics are not generated.
Information accretion (IA) can be calculated as described in
[Wyatt and Radivojac, Bioinformatics, 2013](https://pubmed.ncbi.nlm.nih.gov/23813009/).

```
IDPO:00024  6.32
IDPO:00506  12.04
IDPO:00502  1.34
IDPO:00025  0.56
...
```

## Output files

Output files are generated in the `results` folder. The same files are gerated by both
the command line and the `write_results` function.

* `evaluation_all.tsv` corresponds to the first object returned by the `cafa_eval` function.
* `evaluation_best_< metric >.tsv` corresponds to the second object returned by the `cafa_eval` function. 
A different file for each metric is created.

**Note**: Weighted scores are generated only if the *Information Accretion* file is provided.


## Optional parameters

|   Argument  | Default value | Description                                                                                                                                                                                                                                                |
|:-----------:|------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
|   -out_dir  |   'results' | Output directory (tsv files + log). Either relative to current path or absolute                                                                                                                                                                            |
|     -ia     |            | Information accretion file                                                                                                                                                                                                                                 |
| -no_orphans |  False (flag) | Exclude orphans nodes (e.g. roots) from the calculation                                                                                                                                                                                                    |
|    -norm    |     'cafa'  | Normalization strategy. `cafa` normalize precision by the number of predicted targets and recall by the number of targets in the ground truth. `pred` normalize by the number of  predicted targets. `gt` normalize by the number of ground truth proteins |
|    -prop    |     'max'  | Ancestor propagation strategy. `max` propagate the max score of the traversed subgraph iteratively. `fill` propagate with max until a different score is found                                                                                             |
|   -th_step  |      0.01  | Step size of prediction score thresholds to consider in the range [0, 1). A smaller step, means more calculation                                                                                                                                           |
|  -max_terms |            | Number of terms for protein and namespace to consider in the evaluation. Parsing stops when the target limit for every namespace is reached. The score is not checked, meaning that terms are not sorted before the check, and the check is performed before propagation.                                                                                                                                                                                  |
|   -threads  |       4    | Parallel threads. `0` means use all available CPU threads. Do not use multi thread if you are short in memory                                                                                                                                              |



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

* `fig_< metric >_< name_space >.png` A figure for each namespace in the dataframe and the selected metric. 
The notebook generates one metric at the time, you have to modify the input cell to generate the plots for a different metric

* `fig_< metric >.tsv` A file with the data points for the metric curves. One curve for each method.
```
group	label	ns	tau	cov	wrc	wpr	wf
INGA	INGA_1	biological_process	0.010	0.993	0.557	0.094	0.160
INGA	INGA_1	biological_process	0.020	0.993	0.555	0.094	0.161
INGA	INGA_1	biological_process	0.030	0.993	0.552	0.095	0.162
INGA	INGA_1	biological_process	0.040	0.993	0.551	0.095	0.163
```


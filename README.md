# CAFA-evaluator

This program calculates F-score, weighted F-score and S-score as well as precision-recall and 
remaining uncertainty–misinformation curves as described in the 
[Critical Assessment of protein Function Annotation (CAFA)](https://www.biofunctionprediction.org/cafa/).

CAFA-evaluator is generic and works with any type of ontology. Its implementation is inspired to the 
original Matlab code used in CAFA2 assessment and available at 
[https://github.com/yuxjiang/CAFA2](https://github.com/yuxjiang/CAFA2)

#### CAFA5 challenge
In the Kaggle CAFA5 challenge the software is executed with the following command:

    python3 main.py go-basic.obo predition_dir/ test_terms.tsv -out_dir results/ -ia IA.txt -prop fill -norm cafa -th_step 0.001 -max_terms 500

In the example above the method prediction file should be inside the `prediction_dir/` folder and 
evaluated against the `test_terms.tsv` file containing the ground truth (not available to participants).

## Installation

The software does not require any installation. Simply download the repository and run the **main.py** script. 
The following packages are required:

- numpy
- pandas
- matplotlib (only for generating the plots)

## Algorithm information 

#### Notes
* The word `aspect`, `namespeace` and `sub-ontology` are used interchangeably in the following documentation.
* In order to replicate CAFA results, you can simply adapt the input files. 
  - *No/partial knowledge* can be reproduced by filtering/splitting the **ground truth file** 
  - In order to exclude specific terms from the analyses, 
e.g. generic "binding" terms, you can directly modify the input **ontology file** 

#### Input filtering

Prediction files are filtered considering only those targets included in the ground truth and 
only those terms included in the ontology file. 
If the ground truth contains only annotations from one aspect (e.g. "molecular function"), 
the evaluation is provided only for that aspect.
The `-max_terms` parameter determines the maximum number of terms that will be considered for each target. 
Parsing stops when the target limit for every ontology is reached. The score is not checked, 
meaning that terms are not sorted before the check, and the check is performed before propagation.

The ontology is processed with an internal parser that accepts only the OBO format. 
The following rules are applied: 
  - Obsolete terms are always excluded
  - Only "is_a" and "part_of" relationships are considered
  - Cross-aspect or cross-ontology relationships are always discarded
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



## Assessment

The assessment is provided by running the **main.py** script. While the plots are generated running 
the **plot.ipynb** Jupyter Notebook. To run the assessment:

    python3 main.py ontology.obo predition_dir/ ground_truth_dir/ground_truth.txt -out_dir results/ -ia ia.txt
    
**Mandatory positional arguments**
* *Ontology file* - Ontology in OBO format

* *Prediction folder* - Contain prediction files. Optionally, files can be organized into sub-folders, e.g. one folder per team. 
Sub-folders are processed recursively and the sub-folder name is used as prefix for the method name

* *Ground truth file* - The ground truth file must contain terms in any of the ontology namespace

**Optional arguments (non positional)**

* *-out_dir* - Output folder. Default to current folder
* *-ia* - Information accretion file
* *-prop* - Propagation strategy 
* *-no_orphans* - Exclude orphan terms, e.g. the root(s)
* *-norm* - Normalization strategy

### Input
**Prediction file** - Tab separated file with the target ID, term ID and score columns.

~~~txt
T98230000001    GO:0000159      0.39
T98230000001    GO:0032991      0.39
T98230000001    GO:1902494      0.39
...
~~~

**Ground truth file** - Tab separated file with the target ID and term ID. 
Additional columns are discarded.
~~~
T100900005305   GO:0033234
T100900010085   GO:0006468
T100900010085   GO:0046777
...
~~~

**Information accretion (optional)** - If not provided weighted and S statistics are not generated.
Information accretion (IA) can be calculated as described in
[Wyatt and Radivojac, Bioinformatics, 2013](https://pubmed.ncbi.nlm.nih.gov/23813009/)

```
GO:0000003 3.27
GO:0000035 12.00
GO:0000149 6.82
...
```

### Output


**evaluation_all.tsv** - A table containing the full evaluation, i.e. assessment measures for each threshold. This
file is used as input to generate the plots (see below)
```
filename	ns	tau	cov	pr	rc	f	wpr	wrc	wf	mi	ru	s
INGA_2.cafa	biological_process	0.01000	1.00000	0.07748	0.63499	0.13812	0.05897	0.54060	0.10633	1218.13785	56.96602	1219.46913
INGA_2.cafa	biological_process	0.02000	1.00000	0.07820	0.63330	0.13920	0.05952	0.53867	0.10719	1159.93153	57.13914	1161.33804
...
INGA_1.cafa	cellular_component	0.74000	0.09875	0.58884	0.02458	0.04719	0.58206	0.02551	0.04888	7.93891	297.84363	297.94941
INGA_1.cafa	cellular_component	0.75000	0.05799	0.62162	0.01476	0.02884	0.62162	0.01474	0.02879	7.03334	518.35996	518.40768
...
```

**evaluation_best_< metric >.tsv** - A table containing the best results, best rows from previous file. The metric indicates
based on what metric best rows are selected
```
filename	ns	tau	cov	pr	rc	f	wpr	wrc	wf	mi	ru	s	max_cov
INGA_1.cafa	biological_process	0.50000	0.91253	0.37096	0.34422	0.35709	0.32304	0.28531	0.30300	83.82684	95.50003	127.07161	1.00000
INGA_1.cafa	cellular_component	0.54000	1.00000	0.42442	0.53451	0.47315	0.32350	0.47875	0.38611	26.64830	17.54473	31.90533	1.00000
INGA_1.cafa	molecular_function	0.46000	0.84257	0.52550	0.45361	0.48692	0.50055	0.41148	0.45167	37.42959	30.25438	48.12797	0.99557
INGA_2.cafa	biological_process	0.41000	0.91489	0.36964	0.34685	0.35788	0.32158	0.28766	0.30367	84.40635	94.99940	127.07997	1.00000
INGA_2.cafa	cellular_component	0.52000	1.00000	0.42417	0.53597	0.47356	0.32292	0.48135	0.38653	26.77789	17.38719	31.92757	1.00000
INGA_2.cafa	molecular_function	0.41000	0.92905	0.47084	0.49677	0.48345	0.44646	0.45148	0.44895	43.32693	25.13350	50.08907	0.99778
```


**info.log** - Log file. Information is appended


## Plots

Plots are generated by running all the cells in the **plot.ipynb** Jupyter Notebook after generating the assessment
dataframe **evaluation_all.tsv** (see above). 

### Input

In order to generate the figures you need to manually modify the first cell of the notebook which 
contains information about the input path and a few parameters. Parameters include: 
* the metric (F-score, wighted F-score, S measure)
* the path to the input dataframe (*evaluation_all.tsv*)
* the output folder
* an optional file (*names.tsv*) including information about methods aliases and groups. If provided, the results
are presented selecting only one method per group. Example file:
```
filename	group	label
INGA_1.cafa	INGA	INGA_1
INGA_2.cafa	INGA	INGA_2
```

### Output

**fig_< metric >_< name_space >.png** - A notebook that generates a figure for each namespace in the dataframe and the selected metric. 
The notebook generates one metric at the time, you have to modify the input cell to generate the plots for a different metric

**fig_< metric >.tsv** - A file with the data points for the metric curves. One curve for each method
```
group	label	ns	tau	cov	wrc	wpr	wf
INGA	INGA_1	biological_process	0.010	0.993	0.557	0.094	0.160
INGA	INGA_1	biological_process	0.020	0.993	0.555	0.094	0.161
INGA	INGA_1	biological_process	0.030	0.993	0.552	0.095	0.162
INGA	INGA_1	biological_process	0.040	0.993	0.551	0.095	0.163
```

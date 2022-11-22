# CAFA-evaluator

This program calculates F-score, weighted F-score and S-score as well as precision-recall and 
remaining uncertainty–misinformation curves as described in the 
[Critical Assessment of protein Function Annotation (CAFA)](https://www.biofunctionprediction.org/cafa/).

CAFA-evaluator is generic and works with any type of ontology. It is inspired to the original Matlab code 
used in CAFA2 assessment and available at [https://github.com/yuxjiang/CAFA2](https://github.com/yuxjiang/CAFA2)

In order to replicate CAFA results, you can simply adapt the input files. 
- *No/partial knowledge* can be reproduced by filtering/splitting the **ground truth file** 
- In order to exclude specific terms from the analyses, 
e.g. generic "binding" terms, you can directly modify the input **ontology file** 


## Algorithm information 

#### Input filtering

Prediction files are filtered considering only those targets included in the ground truth and those terms included in 
the ontology file. If the ground truth contains only annotations from one aspect (e.g. molecular function), 
the evaluation is provided only for that aspect

The ontology is processed with an internal parser. Only the OBO format is allowed. The following is also 
  - Obsolete terms are always excluded
  - Only "is_a" and "part_of" relationships are considered
  - Cross-aspect or cross-ontology relationships are always discarded
  - Alternative term IDs are automatically mapped to the main ID

When information accretion is provided, terms which are not available in the file are removed from the ontology


#### Terms propagation

Both the predictions and the ground truth annotations are always propagated up to the ontology root(s). By
default, prediction scores are propagated without overwriting the scores assigned to the parents. Alternatively,
an option let you to propagate considering always the max

#### Memory usage

- The algorithm stores in memory a Numpy boolean *N x M* array 
(N = number of ground truth targets; M = ontology terms of a single aspect)
for each aspect in the ground truth file.

- An array of the same size (rows &le; N), but containing floats (the prediction scores) instead of booleans, 
is stored for each prediction file. Prediction files are processed one by one and the matrix gets reassigned.


#### Required packages

- numpy
- pandas
- matplotlib (only for the plots)

## Assessment

The assessment is provided by running the **main.py** script. While the plots are generated running 
the **plot.ipynb** Jupyter Notebook. To run the assessment:

    python3 main.py ontology.obo predition_dir/ grounf_truth_dir/ground_truth.txt -out_dir results/ -ia ia.txt
    
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


**hmean_< metric >.tsv** - The harmonic mean across namespace. Example output file *eval_f.tsv*
```
filename	rc	pr	cov	max_cov	f
INGA_1.cafa	0.443	0.419	0.946	0.999	0.431
INGA_2.cafa	0.444	0.417	0.947	0.999	0.430
...
```

**df_all.tsv** - A dataframe containing the full evaluation, i.e. assessment measures for each threshold. This
file is used as input to generate the plots (see below)
```
filename	ns	tau	cov	pr	rc	f	wpr	wrc	wf	mi	ru	s
INGA_2.cafa	biological_process	0.01	1.0	0.077484756570578	0.6349876091296363	0.13811584192570855	0.058965195042090884	0.5405965440442712	0.10633227099255976	1218.1378520809258	56.96601624021818	1219.4691278087412
INGA_2.cafa	biological_process	0.02	1.0	0.07819651116457993	0.6333035302891409	0.13920484523272886	0.05951535470382788	0.5386713751192704	0.10718799451971453	1159.931532153873	57.13913838968864	1161.3380388245048
INGA_2.cafa	biological_process	0.03	1.0	0.08052587562894702	0.625829378101304	0.14269153772068444	0.061112601760945184	0.5322187141176382	0.10963611545576801	1113.1721207892651	57.875847431198515	1114.67564081142
...
```


**info.log** - Log file. Information is appended


## Plots

Plots are generated by running all the cells in the **plot.ipynb** Jupyter Notebook after generating the assessment
dataframe **df_all.tsv** (see above). 

### Input

In order to generate the figures you need to manually modify the first cell of the notebook which 
contains information about the input and a few parameters. Parameters include: 
* the metric (F-score, wighted F-score, S measure)
* the path to the input dataframe (*df_all.tsv*)
* the output folder
* an optional file (*names.tsv*) including information about method alias and group. If provided, the results
are presented selecting only one method per group. For example:
```
filename	group	label
INGA_1.cafa	INGA	INGA_1
INGA_2.cafa	INGA	INGA_2
```

### Output

**fig_< metric >_< name_space >.png** - The notebook a figure for each namespace in the dataframe and selected metric. 
The notebook generates one metric at the time, you have to modify the input cell to generate the plots for a different metric

**eval_< metric >.tsv** - A file with the data points for the metric curves. One curve for each method
```
group	label	ns	tau	cov	wrc	wpr	wf
INGA	INGA_1	biological_process	0.010	0.993	0.557	0.094	0.160
INGA	INGA_1	biological_process	0.020	0.993	0.555	0.094	0.161
INGA	INGA_1	biological_process	0.030	0.993	0.552	0.095	0.162
INGA	INGA_1	biological_process	0.040	0.993	0.551	0.095	0.163
```

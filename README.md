# CAFA-evaluator

Calculate F-score and S-score as well as precision-recall and remaining uncertainty–misinformation curves  
as provided in the 
[Critical Assessment of protein Function Annotation (CAFA)](https://www.biofunctionprediction.org/cafa/).
CAFA-evaluator is generic and works with any type of ontology.



### Algorithm information 

Important information about the algorithm.

- Prediction files are filtered considering only targets and ontology terms which are included 
in the ground truth and ontology files.

- Only aspects which are included in the ground truth file are considered. For example, if the ground truth contains 
only molecular function annotations, the output will only cover molecular function analyses, 
even if predictions for other aspects are provided.

- Both the predictions and the ground truth annotations are always propagated up to the ontology root(s).

- The algorithm stores in memory a Numpy boolean *N x M* array (N = number of ground truth targets; M = ontology terms of a single aspect)
for each aspect in the ground truth file.

- An array of the same size (rows &le; N), but containing floats (the prediction scores) instead of booleans, 
is stored for each prediction file. 
Prediction files are processed one by one and the matrix gets reassigned.

- The ontology is processed with an internal parser.
  - Only the OBO format is allowed.
  - Obsolete terms are always excluded.
  - Cross-aspect or cross-ontology relationships are always discarded.
  - Only "is_a" and "part_of" relationships are considered.
  - Alternative term IDs are automatically mapped to the main ID.

- Information accretion, when provided, is automatically set to zero for "nan" and "inf" values.

In order to replicate CAFA results, you can simply adapt the input files. 
- By filtering/splitting the input **ground truth file** you can replicate the "no knowledge" and "partial knowledge" benchmarks. 
- In order to exclude specific terms from the analyses, e.g. generic "binding" terms, you can directly modify the input **ontology file**. 


## Required packages

- numpy
- matplotlib
- pandas

## Execution

To run CAFA-evaluator simply call **main.py**.

    python3 main.py pontology.obo pred_dir/ gt_dir/gt.txt -out_dir results/ -ia ia.txt
    

**Mandatory arguments** (positional)
* *Ontology file* - Ontology in OBO format

* *Prediction folder* - Contain prediction files. Optionally, files can be organized into sub-folders, e.g. one folder per team. 
Sub-folders are processed recursively and the sub-folder name is used as prefix for the method name

* *Ground truth file* - The ground truth file must contain terms in any of the ontology namespace

**Optional** (non positional)

* *-out_dir* - Output folder. Default to current folder
* *-names* - File with information about methods (filename, group, label, is_baseline). Used to select only one method per team, 
and for the labels in the plots. If not provided all methods are considered and labels are the file names 
* *-ia* - Information accretion file

## Input format
**Prediction**

Tab separated file with the target ID, term ID and score columns.

~~~txt
T98230000001    GO:0000159      0.39
T98230000001    GO:0032991      0.39
T98230000001    GO:1902494      0.39
...
~~~

**Ground truth**

Tab separated file with the target ID and term ID. 
Additional columns are discarded.
~~~
T100900005305   GO:0033234
T100900010085   GO:0006468
T100900010085   GO:0046777
...
~~~

**Names (optional)**

Tab separated file with header.
This is used to filter the best method for each team and 
set labels in the output figures. 
If not provided, the program will use all predictions and 
the labels will be file names.

```
name   group   label   is_baseline
naive_pred_1.txt    naive   Naive   True
blast_pred_1.txt    blast   Blast   True
team_8_pred_1.txt   team_8  team8_m1    False
team_4_pred_1.txt   team_4  team4_m1    False
```

**Information accretion (optional)**

If not provided the corresponding plot files are not generated.
Information accretion (IA) can be calculated as described in 
[Wyatt and Radivojac, Bioinformatics, 2013](https://pubmed.ncbi.nlm.nih.gov/23813009/).


```
GO:0000003 3.27
GO:0000035 12.00
GO:0000149 6.82
GO:0000322 4.83
GO:0000323 3.11
GO:0000325 5.70
GO:0000407 8.15
GO:0000408 9.67

```


### Output

#### fig_< metric >_< name_space >.png
The metric can be: 
  * *f*, Precision / recall curves
  * *wf*, Weighted precision / recall curves. This file is generated only when IA is provided.
  * *s*, Misinformation / remaining uncertainty curves. This file is generated only when IA is provided.

#### eval_< metric >.tsv - The corresponding curve points
Example output file *eval_f.tsv*:
```
group	label	ns	tau	cov	rc	pr	f
INGA	INGA_2	biological_process	0.010	0.993	0.641	0.117	0.197
INGA	INGA_2	biological_process	0.020	0.993	0.640	0.117	0.198
INGA	INGA_2	biological_process	0.030	0.993	0.637	0.118	0.199
INGA	INGA_2	biological_process	0.040	0.993	0.635	0.119	0.201
INGA	INGA_2	biological_process	0.050	0.993	0.633	0.121	0.203
INGA	INGA_2	biological_process	0.060	0.993	0.629	0.124	0.208
INGA	INGA_2	biological_process	0.070	0.993	0.629	0.124	0.208
INGA	INGA_2	biological_process	0.080	0.993	0.626	0.126	0.210
```

#### hmean_< metric >.tsv - The harmonic mean across namespace
Example output file *eval_f.tsv*:
```
group	cov	max_cov	f
INGA	0.943	0.992	0.470
```

#### info.log
Log file. Information is appended
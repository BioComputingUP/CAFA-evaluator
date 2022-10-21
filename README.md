# CAFA-evaluator

Calculate F-score and S-score as well as precision-recall and remaining uncertainty–misinformation curves  
as provided in the 
[Critical Assessment of protein Function Annotation (CAFA)](https://www.biofunctionprediction.org/cafa/)

CAFA-evaluator is generic and works with any ontology. For example, it can be used to evaluate 
intrinsic disorder function predictions against [DisProt](https://disprot.org/) annotations.



## Required packages

- numpy
- logging
- matplotlib

## HowTo

To run CAFA-evaluator simply call **main.py**.

    python3 main.py pontology.obo pred_dir/ gt_dir/ -out_dir results/ -methods methods.txt -ia ia.txt
    

**Mandatory arguments** (positional)
* *Ontology file* - Ontology in OBO format

* *Prediction folder* - Contain prediction files. Optionally, files can be organized into sub-folders, e.g. one folder per team. 
Sub-folders are processed recursively and the sub-folder name is used as prefix for the method name

* *Ground truth folder* - Contain ground truth files, one for each namespace (sub-ontology). 
The ontology namespace has to be hardcoded in the filename, e.g. gt_*bpo*_.txt or *disorderfunction*.tsv

**Optional** (non positional)

* *-out_dir* - Output folder. Otherwise it will use the current folder
* *-filter_dir* - Benchmark folder with benchmark files, one for each namespace (sub-ontology). 
The evaluation is performed considering only the targets listed in these files. If not provided all target in
the ground truth are considered.
The ontology namespace has to be hardcoded in the filename, e.g. bm_*bpo*_.txt or *disorderfunction*.tsv
* *-no_roots* - Exclude terms without is_a relationships (roots). Use with caution if the ontology has been filtered
* *-names* - File with information about methods (filename, group, label, is_baseline). Used to select only one method per team, 
and for the labels in the plots. If not provided all methods are considered and labels are the file names 
* *-ia* - Dictionary with information accretion




### Namespaces (sub-ontologies)

At the moment in order to recognize the namespaces you should 
modify this dictionary which is hardcoded in *main.py*.

The key is the string matched in the ground truth (and benchmark) files, 
while values correspond to the namespace field in the OBO file. 

~~~python3
namespaces = {
    "bpo": "biological_process", 
    "cco": "cellular_component", 
    "mfo": "molecular_function",
    "disorderfunction": "Disorder function", 
    "interactionpartner": "Interaction partner", 
    "structuralstate": "Structural state", 
    "structuraltransition": "Structural transition"
}
~~~

### Output

* prrc_<name_space>.png - Precision-recall plots, one for each ontology
* prrc_<name_space>.tsv - Precision-recall points, one for each ontology
* miru_<name_space>.png - Remaining uncertainty–misinformation plots, one for each ontology
* miru_<name_space>.tsv - Remaining uncertainty–misinformation points, one for each ontology
* info.log - Log file. Information is appended

prrc_<name_space>.tsv
```
ns      group   method  tau     pr      rc      f       cov_f   mi      ru      s       cov_s   name    label   is_baseline
Disorder_function       team_7  team_7_pred_1.txt       0.990   0.556   0.018   0.035   0.018   0.101   51.885  51.885  0.018   team_7_pred_1.txt       team7_m1        False
Disorder_function       team_7  team_7_pred_1.txt       0.980   0.556   0.018   0.035   0.018   0.101   51.885  51.885  0.018   team_7_pred_1.txt       team7_m1        False
Disorder_function       team_7  team_7_pred_1.txt       0.970   0.556   0.018   0.035   0.018   0.101   51.885  51.885  0.018   team_7_pred_1.txt       team7_m1        False
```

miru_<name_space>.tsv
```
ns      group   method  tau     pr      rc      f       cov_f   mi      ru      s       cov_s   name    label   is_baseline
Disorder_function       team_6  team_6_pred_1.txt       0.640   1.000   0.003   0.006   0.003   0.000   0.000   0.000   0.003   team_6_pred_1.txt       team6_m1        False
Disorder_function       team_6  team_6_pred_1.txt       0.630   1.000   0.003   0.006   0.003   0.000   0.000   0.000   0.003   team_6_pred_1.txt       team6_m1        False
Disorder_function       team_6  team_6_pred_1.txt       0.620   1.000   0.003   0.006   0.003   0.000   0.000   0.000   0.003   team_6_pred_1.txt       team6_m1        False
```

### Prediction file format

Target_ID   Term_ID Score

~~~txt
T98230000001    GO:0000159      0.39
T98230000001    GO:0032991      0.39
T98230000001    GO:1902494      0.39
...
~~~

### Ground truth file format

A file containing the Target ID and Term ID as below.
~~~
T100900005305   GO:0033234
T100900010085   GO:0006468
T100900010085   GO:0046777
...
~~~


### Benchmark file format

Files, one for each ontology containing the list of protein IDs.
~~~
T446890002807
T446890001819
T446890003376
...
~~~

### Names file format (tab separated)

It must be a TAB separated file with header.
This is used to consider the best method for each team and set labels in the output figures. 
If this file is not provided, the program will use all predictions and 
the labels will be the file names.

```
name   group   label   is_baseline
naive_pred_1.txt    naive   Naive   True
blast_pred_1.txt    blast   Blast   True
team_8_pred_1.txt   team_8  team8_m1    False
team_4_pred_1.txt   team_4  team4_m1    False
```

### Information accretion file format

Information accretion (IA) can be calculated as described in 
[Wyatt and Radivojac, Bioinformatics, 2013](https://pubmed.ncbi.nlm.nih.gov/23813009/).
It must be a TAB separated file with header. Only the *term* and *ia* columns are used.


```
term    co_occuring     co_occuring_parents     p_cond  ia
DO:00000        484     484     1.0     -0.0
DO:00001        141     484     0.29132231404958675     1.7793118848758012
DO:00002        106     141     0.75177304964539        0.4116308978355944
DO:00003        5       141     0.03546099290780142     4.817623257511431
DO:00005        2       141     0.014184397163120567    6.139551352398794
DO:00007        2       141     0.014184397163120567    6.139551352398794
DO:00008        154     484     0.3181818181818182      1.6520766965796931
```

# Change Log
All notable changes to this project will be documented in this file.
 
## [Unreleased] - 2024-01-23
  
We changed the way alternative identifiers in ontology files are considered.
Now alternative identifiers are recognized in both the ground truth and prediction files
and mapped to the "canonical" term.
 
### Added
- CHANGELOG.md, this file!
 
### Changed
- graph.py, changed the Graph class.
- parser.py, changed the ground truth and prediction parsers in order 
to replace alternative identifiers with canonical terms.
- plot.ipynb, cleaned up.

 
## [1.0.0] - 2023-08-04
 
First release after CAFA5 challenge closed. The version used in Kaggle is
provided under the 'Kaggle' branch. The 'main' branch includes instead a
packed version of the code. While usage and performance has been improved 
the calculation is exactly the same.
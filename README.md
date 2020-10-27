# Assay targeted hERG

## Dependency
- Python 3.7.x
    - library RDKIT on docker
    - numpy, scipy

- R > 3.6

- From personal scripts
	- QSAR R scripts from github => https://github.com/ABorrel/QSAR-QSPR
	- Toolbox R => https://github.com/ABorrel/R_toolbox
	- CompDesc library (pip install -i https://test.pypi.org/simple/ CompDesc)

## To do
- ~~Add function to merge QSAR model from undersampling~~
- ~~SOM add descriptors significance~~ 
- Build external dataset on chembl for active chemicals on herg
- merge AD scripts with the QSAR-QSPR scripts


## Update
- 08-6-20: Unsupervised analysis 
- 10-6-20: Add QSAR-QSPR model builder
- 16-6-20: merge models from different machine learning
- 1-7-20: Fix error in R runner 
- 6-7-20: correct SOM when apply
- 15-7-20: Add AD computation
- 7-8-20: fix AD for figure and push clustering
- 31-8-20: add model regression
- 05-10-20: change hyperparameter for regression modeling

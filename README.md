SAGA
====

Software for the Analysis of Genetic Architecture

to cite this package please use:

Blackmon, H. and J. P. Demuth. 2015. SAGA: An R package for an information-theoretic approach to estimating the composite genetic effects contributing to variation among generation means: moving beyond the joint-scaling test for line cross analysis.  [![DOI:10.5281/zenodo.32513](https://zenodo.org/badge/18367/coleoguy/SAGA.svg)](https://zenodo.org/badge/latestdoi/18367/coleoguy/SAGA)




The package SAGA is a collection functions to ease and hopefully improve the quality of line cross analysis of genetic architecture.  The overall goal is to allow for an easy and straightforward implementation of model averaged analysis using AICc.

If you want the very latest version of SAGA then you can install my developmental version from GitHub.

If you don't have the package devtools installed yet you will need to do that first: 

`install.packages("devtools")`

`library(devtools)`

`install_github('coleoguy/SAGA')`

`library(SAGA)`


This package also includes a vignette to build the vignette during instillation use the following code:

`install_github("coleoguy/SAGA", build_vignettes = TRUE)`

`library(SAGA)`

If you are using R studio the best way to get to the vignette is to click on the package in the packages tab in the lower
right window.  Then there will ba a link to package descriptions and vignettes. If your just using a bare terminal you can use this command:

`vignette("model-averaged-analysis", package='SAGA')`

if you have questions or problems please let me know [coleoguy@gmail.com](mailto:coleoguy@gmail.com).

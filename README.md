# Bayes Optimal Informer SEts

This is Github repo for BOISE package as well as codes for accelerated BOISE and missing data approach. 

## Main functions in BOISE

* dpmm_beta(): DPMM clustering for binary interaction data. 

* pel1_beta(): Compute PEL$_1$ loss for a given informer set.

* Boise(): Bayes optimal informer set selection method. Two stage decision procedure with adaptive selection. 

* Boise_Aug(): Augment an existing informer set with BOISE method same as above. 

* Evaluate(): Evaluate an informer set on the test set. The metric could be one of NEF, ROCAUC, MCC, F1 score.

## For BOISE R package installation

```R
install.packages("devtools")
devtools::install_github("wiscstatman/esdd/BOISE")
```

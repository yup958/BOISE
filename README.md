# DPMM-for-Rasch
Apply Dirichlet Process Mixture Model on Rasch Model Clustering in R.

Our data is a binary matrix with 200 targets (rows) by 366 compounds (columns). Our goal is to learn a clustering structure on targets (rows).

Our model setting is that whole matrix comes from a rasch model, with targets in the same clusters share same parameters (alphas and betas). 

This is an initial step in Bayes Optimal Informer SEt (BOISE).

Technically speaking, we follow the Algorithm 8 from paper *Markov Chain Sampling Methods for Dirichlet Process Mixture Model, R.M.Neal*, with non-conjugate prior and auxillary variable method.

There are 9 functions: 

* prob() for calculating density;

* parasample() for auxillary variables sampling;

* Initialization() for giving initial cluster assignment from CRP;

* Rearrange() for rearrange clusters based on auxillary variable method;

* Rasch_Matrix.stan and Rasch_Vector.stan are stan file to get a posterior estimation of parameters after clustering;

* par_update() for updating parameters after rearrangement;

* label_update() for delete empty clusters after rearrangement;

* dpmm() to make everything work together. There is a code test in __dpmm.R__, with 10 iterations, concentration parameter a = 20, auxillary number 10. It will take approximate 3 hours...

This is a very naive version...

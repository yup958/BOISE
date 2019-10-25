# BOISE with Binomial Model and Beta Conjugate Prior

Codes within this directory are based on Bernoulli-Beta model, with conjugate prior and apply Algorithm 3 in Neal's 2000 paper.

## BOISE

There are 2 functions directly related to BOISE estimation:

* pel2_beta() for calculating posterior probability p(xA|x0, Cl), PEL2 value and top set for a given clustering assignment cl, given intermediate data xA.

* pel1_beta() for calculating pel1 value for given informer set A and priors, for a bunch of DPMM clustering samples cl_sample.

## DPMM for Conjugate Beta Prior

Apply Dirichlet Process Mixture Model on Binomial-Beta Model Clustering in R.

Our data is a binary matrix with 200 targets (rows) by 366 compounds (columns). Our goal is to sample a bunch of clustering assignments on targets (rows).

Our model setting is that whole matrix comes from a Binomial model, with targets in the same clusters share same success probability (\phi_{kj}). 

This is an initial step in Bayes Optimal Informer SEt (BOISE).

Technically speaking, we follow the Algorithm 3 from paper *Markov Chain Sampling Methods for Dirichlet Process Mixture Model, R.M.Neal*, with conjugate prior and collapsing the parameter \phi.

There are 3 functions: 

* Initial_beta() for giving initial cluster assignment from CRP;

* Update_beta() for updating cluster assignments once based on Algorithm 3;

* dpmm_beta() for sampling a bunch of cluster assignments, after the initial warm-up process;

* Analysis is a R script for test this model. Could be update

Typically, we use hyperparameter alpha and beta to maintain the overall average hit rate within each colomn (compound).
This is a 2nd version of BOISE model...


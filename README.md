# DPMM-for-Rasch
Apply Dirichlet Process Mixture Model on Rasch Model Clustering.
Our data is a binary matrix with 200 targets (rows) by 366 compounds (columns). Our goal is to learn a clustering structure on targets (rows).
Our model setting is that whole matrix comes from a rasch model, with targets in the same clusters share same parameters (alphas and betas). 
This is an initial step in Bayes Optimal Informer SEt (BOISE).

Technically speaking, we follow the Algorithm 8 from paper *Markov Chain Sampling Methods for Dirichlet Process Mixture Model, R.M.Neal*, with non-conjugate prior and auxillary variable method.
This is a very naive version...

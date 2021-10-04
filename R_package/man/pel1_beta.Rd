\name{pel1_beta}
\alias{pel1_beta}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  function to compute PEL-1 loss of an informer set
}
\description{
Input a sample of clustering assignments cl_sample, sample size iter, informer set A, list P, with P[[i]] a matrix summarizing the ith cluster assginment. top size nT, prior alpha and beta, data x0, divergence m0
Output the PEL1 of A
}
\usage{
pel1_beta(cl_sample, P, sample_size, interm_size, A, nT = 10, alpha, beta, x0, m0 = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cl_sample}{
  An augmented list of clustering assignments as well as the sampled outcomes on the new target i*. cl_sample now is {KK, NN, CC, XX}, where XX contains "size" samples of possible outcomes on i*.
}
  \item{P}{
  An intermediate matrix P that is computed by clust_sum function. It summarizes the cluster assignment information of given cl_sample.
}
  \item{sample_size}{
Sample size of DPMM clustering
}
  \item{interm_size}{
Sample size of possible intermediate outcomes xA
}
  \item{A}{
Informer set
}
  \item{nT}{
Length of top set, used to define loss function in BOISE.
}
  \item{alpha}{
Prior hyperparameters alpha0
}
  \item{beta}{
Prior hyperparameters beta0
}
  \item{x0}{
Initial bioactivity data.
}
  \item{m0}{
Chinese restaurant prior mass m0
}
}
\details{
alpha, beta, m0 should keep the same as in \code{dpmm_beta}. nT is the top set size that defines loss function. P should come from function \code{clust_sum}. 
}
\value{
Output is the posterior expected loss (PEL1) of certain informer set A. Defined in BOISE paper.
}
\references{
Yu, P., Ericksen, S.S., Gitter, A. and Newton, M.A. (2020) \emph{Bayes Optimal Informer Sets for Drug Discovery}. University of Wisconsin, Technical Report, 2020.
}
\author{
Peng Yu, Michael A. Newton
}
\note{
This is an internal function used in function \code{Boise}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\url{https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006813}
}
\examples{
#suppose you have the data and cl_sample as in dpmm_beta
P = clust_sum(cl_sample, x0, sample_size = length(cl_sample$KK), alpha, beta)
sample_size = length(cl_sample$KK)
interm_size = 1000
m = nrow(x0)
n = ncol(x0)
cl_sample$XX = rep(0, sample_size * interm_size * n)
dim(cl_sample$XX) = c(sample_size, interm_size, n)
for (j in 1:sample_size) {
    toy_samples = rbinom(interm_size * n, 1, 0.5)
    dim(toy_samples) = c(interm_size, n)
    cl_sample$XX[j, , ] = toy_samples
}
pel1 = pel1_beta(cl_sample, P,sample_size, interm_size, A = c(1,2), nT=10, alpha, beta, x0, m0 = 15)
}
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

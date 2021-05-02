\name{pel2_beta}
\alias{pel2_beta}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  function to calculate PEL-2 loss based on intermediate data xA
}
\description{
Input: Matrix P summarize the cluster assignment information (For computation)
        Informer set A, Data matrix x0, Intermediate data xA, Top size nT, 
        prior alpha, beta, divergence m0
Output: PEL2 value (given x0,xA,Cl) [posterior expected theta E(theta|C,x0,xA)]
}
\usage{
pel2_beta(P, x0, xA, A, nT = 10, alpha, beta, m0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{P}{
  An intermediate matrix P that is computed by clust_sum function. It summarizes the cluster assignment information of given cl_sample.
}
  \item{x0}{
  Initial bioactivity data.
}
  \item{xA}{
  Intermediate data on the selected informer set
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
  \item{m0}{
  Chinese restaurant prior mass m0
}
}
\details{
Here P is not the list of matrices returned by function \code{clust_sum}. It stands for one such matrix in the list.
}
\value{
  PEL2 value (given x0,xA,Cl) [posterior expected theta E(theta|C,x0,xA)]
}
\references{
Yu, P., Ericksen, S.S., Gitter, A. and Newton, M.A. (2020) \emph{Bayes Optimal Informer Sets for Drug Discovery}. University of Wisconsin, Technical Report, 2020.
}
\author{
Peng Yu, Michale A. Newton
}
\note{
This is an internal function used in function \code{pel1_beta}, also for function \code{Boise}.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\url{https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006813}
}
\examples{
#suppose you have the dat and cl_sample as in dpmm_beta
P = clust_sum(cl_sample, x0, sample_size = length(cl_sample$KK), alpha, beta)
xA = c(0,1,0,1,0)
A = c(1,2,3,4,5)
nT = 10
post_theta = pel2_beta(P[[1]], x0, xA, A, nT = 10, alpha, beta, m0 = 15)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

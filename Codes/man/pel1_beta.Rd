\name{pel1_beta}
\alias{pel1_beta}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  function to compute PEL-1 loss of an informer set
}
\description{
Input a sample of clustering assignments cl_sample, sample size iter, informer set A and nA, list P, with P[[i]] a matrix summarizing the ith cluster assginment. top size nT, prior a and b, data x0, divergence alpha
Output the PEL1 of A
}
\usage{
pel1_beta(cl_sample, P, iter, size, A, nA, nT = 10, a, b, x0, alpha = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cl_sample}{
  An augmented list of clustering assignments as well as the sampled outcomes on the new target i*. cl_sample now is {KK, NN, CC, XX}, where XX contains "size" samples of possible outcomes on i*.
}
  \item{P}{
  An intermediate matrix P that is computed by clust_sum function. It summarizes the cluster assignment information of given cl_sample.
}
  \item{iter}{
Sample size of DPMM clustering
}
  \item{size}{
Sample size of possible intermediate outcomes xA
}
  \item{A}{
Informer set
}
  \item{nA}{
Length of informer set
}
  \item{nT}{
Length of top set, used to define loss function in BOISE.
}
  \item{a}{
Prior hyperparameters alpha0
}
  \item{b}{
Prior hyperparameters beta0
}
  \item{x0}{
Original interaction data.
}
  \item{alpha}{
Chinese restaurant prior mass m0
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
#suppose you have the dat and cl_sample as in dpmm_beta
P = clust_sum(cl_sample, dat, iter = length(cl_sample$KK), a, b)
size = 1000
n = nrow(dat)
m = ncol(dat)
cl_sample$XX = rep(0, iter * size * m)
dim(cl_sample$XX) = c(iter, size, m)
for (j in 1:iter) {
    toy_samples = rbinom(size*m, 1, 0.5)
    dim(toy_samples) = c(size, m)
    cl_sample$XX[j, , ] = toy_samples
}
pel1 = pel1_beta(cl_sample, P,iter, size, A = c(1,2), nA=length(A), nT=10, a, b, x0 = dat, alpha = 15)
}
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

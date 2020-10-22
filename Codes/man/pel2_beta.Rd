\name{pel2_beta}
\alias{pel2_beta}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  function to calculate PEL-2 loss based on intermediate data xA
}
\description{
Input: Matrix P summarize the cluster assignment information (For computation)
        Informer set A, Data matrix x0, Intermediate data xA, 
        Informer size nA, Top size nT, 
        prior a,b, divergence alpha
Output: PEL2 value (given x0,xA,Cl) [posterior expected theta E(theta|C,x0,xA)]
}
\usage{
pel2_beta(P, x0, xA, nA, A, nT = 10, a, b, alpha)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{P}{
  An intermediate matrix P that is computed by clust_sum function. It summarizes the cluster assignment information of given cl_sample.
}
  \item{x0}{
  Original interaction data.
}
  \item{xA}{
  Intermediate data on the selected informer set
}
  \item{nA}{
  Length of informer set
}
  \item{A}{
  Informer set
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
  \item{alpha}{
  Chinese restaurant prior mass m0
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  PEL2 value (given x0,xA,Cl) [posterior expected theta E(theta|C,x0,xA)]
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
xA = c(0,1,0,1,0)
A = c(1,2,3,4,5)
nT = 10
post_theta = pel2_beta(P[[1]], dat, xA, nA=length(A), A, nT = 10, a, b, alpha = 15)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

\name{dpmm_beta}
\alias{dpmm_beta}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  function to do MCMC sampling for conjugate prior (Beta-binomial) clustering assignments
}
\description{
Output: A list cl_sample, including clustering assignments of size "iter"
}
\usage{
dpmm_beta(a, b, x0, warm = 200, iter = 10, step = 5, alpha = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{a}{
  Prior hyperparameters alpha0
}
  \item{b}{
  Prior hypereparameters beta0
}
  \item{x0}{
  Original interaction data x0
}
  \item{warm}{
  Burn-in parameter. How many samples are dropped in the beginning.
}
  \item{iter}{
  Sample size of DPMM
}
  \item{step}{
  Thinning parameters. Gap between consecutive samples.
}
  \item{alpha}{
  Chinese restaurant prior mass m0.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
cl_sample, is a list of clustering assignment samples. cl_sample = {KK, NN, CC}, where KK is an array  of size "iter" with each element as the number of clusters; NN is a "iter" * "row(x0)" matrix with each row as the cluster assignment of the targets; CC is summary of every cluster in each assignment. 
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
data(pkis1)
## Use 2sd criteria to create binary matrix
dat <- t(apply(pkis1, 1, function(x){
  thres = mean(x) + 2 * sd(x)
  return(as.numeric(x>thres))
}))
a = rep(mean(dat,na.rm = T),ncol(dat))
b = 1-a
cl_sample = dpmm_beta(a,b,x0 = dat, warm=500, iter=10, step=10, alpha=15)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

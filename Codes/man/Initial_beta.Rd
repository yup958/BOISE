\name{Initial_beta}
\alias{Initial_beta}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  intermediate function to initial Chinese Restaurant Process for conjugate prior DPMM
}
\description{
Output a list of {K, N, C}, where K = number of clusters, N = #of targets in each cluster, C = cluster assignment.
}
\usage{
Initial_beta(dat, a, b, alpha = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
  The original interaction data x0
}
  \item{a}{
  Beta prior hyperparameters alpha0 
}
  \item{b}{
  Beta prior hyperparameters beta0
}
  \item{alpha}{
  Chinese restaurant process prior mass m0.
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
data(pkis1)
## Use 2sd criteria to create binary matrix
dat <- t(apply(pkis1, 1, function(x){
  thres = mean(x) + 2 * sd(x)
  return(as.numeric(x>thres))
}))
rm(pkis1)
a = rep(mean(dat),ncol(dat))
b = 1-a
cl = Initial_beta(dat,a,b,10)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

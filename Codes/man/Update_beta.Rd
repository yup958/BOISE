\name{Update_beta}
\alias{Update_beta}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Dirichlet process updating for one iteration (based on Neal's 2000 Paper, Algorithm 3)
}
\description{
Update the clustering assignment list cl = {K, N, C}
}
\usage{
Update_beta(cl, x0, alpha, beta, m0 = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cl}{
  Old clustering assignment list of {K, N, C}
}
  \item{x0}{
  Initial bioactivity data x0
}
  \item{alpha}{
  Hyperparameters alpha0
}
  \item{beta}{
  Hyperparameters beta0
}
  \item{m0}{
  Chinese restaurant prior mass m0
}
}
\details{
alpha, beta and m0 are hyperparameters determined by user. x0 is the original data used throughout the process. 
}
\value{
An updated clustering assignment list cl = {K, N, C}, where K is the number of clusters, N is number of targets in each cluster, and C is the assignment of each target.
}
\references{
Yu, P., Ericksen, S.S., Gitter, A. and Newton, M.A. (2020) \emph{Bayes Optimal Informer Sets for Drug Discovery}. University of Wisconsin, Technical Report, 2020.}
\author{
Peng Yu, Michael A. Newton
}
\note{
This is an internal function for \code{dpmm_beta}.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\url{https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006813}
}
\examples{
data(pkis1)
## Use 2sd criteria to create binary matrix
x0 <- t(apply(pkis1, 1, function(x){
  thres = mean(x) + 2 * sd(x)
  return(as.numeric(x>thres))
}))
rm(pkis1)
alpha = rep(mean(x0,na.rm = T),ncol(x0))
beta = 1-alpha
m0 = 10
cl = Initial_beta(x0,m0)
cl = Update_beta(cl, x0, alpha, beta, m0)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

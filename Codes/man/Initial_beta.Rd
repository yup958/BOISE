\name{Initial_beta}
\alias{Initial_beta}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Internal function to initial Chinese Restaurant Process for conjugate prior Dirichlet process mixture model.
}
\description{
Output a list of number of clusters, number of targets in each cluster, and cluster assignment.
}
\usage{
Initial_beta(x0, m0 = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x0}{
  The initial bioactivity data x0
}
  \item{m0}{
  Chinese restaurant process prior mass m0.
}
}
\details{
m0 is a hyperparameters determined by users. 
}
\value{
cl_sample, a list of {K, N, C}, where K = number of clusters, N = #of targets in each cluster, C = cluster assignment.
}
\references{
Yu, P., Ericksen, S.S., Gitter, A. and Newton, M.A. (2020) \emph{Bayes Optimal Informer Sets for Drug Discovery}. University of Wisconsin, Technical Report, 2020.
}
\author{
Peng Yu, Michael A. Newton
}
\note{
This is just an internal function for \code{dpmm_beta}
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
cl = Initial_beta(x0,m0=10)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

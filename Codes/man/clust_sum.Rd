\name{clust_sum}
\alias{clust_sum}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  An intermediate function that summarizes the cl_sample to P.
}
\description{
  Summary the clustering information for each sample from MCMC. Including number of clusters, number of targets in each cluster, and the actives of each compound within each cluster.
}
\usage{
clust_sum(cl_sample, x0, sample_size, alpha, beta)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cl_sample}{
  A list that contains all the clustering assignment information
}
  \item{x0}{
  Initial bioactivity data.
}
  \item{sample_size}{
  Sample size of DPMM clustering samples
}
  \item{alpha}{
  Prior hyperparameters alpha0
}
  \item{beta}{
  Prior hyperparameters beta0
}
}
\details{
x0 is the original data used throughout the process. alpha, beta and sample_size are the same as the corresponding parameters in function \code{dpmm_beta}.
}
\value{
P, a list of matrices, with P[[i]] as a matrix contain all information of ith cluster assignment. P[[i]] is K by (n+1) matrix. Each row for one cluster, last column as number of targets in the cluster.
}
\references{
Yu, P., Ericksen, S.S., Gitter, A. and Newton, M.A. (2020) \emph{Bayes Optimal Informer Sets for Drug Discovery}. University of Wisconsin, Technical Report, 2020.
}
\author{
Peng Yu, Michael A. Newton
}
\note{
Pleas test after you have the data and cl_sample created after dpmm_beta.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\url{https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006813}
}
\examples{
#suppose you first have the dat and cl_sample as in dpmm_beta
P = clust_sum(cl_sample, x0, sample_size = length(cl_sample$KK), alpha, beta)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

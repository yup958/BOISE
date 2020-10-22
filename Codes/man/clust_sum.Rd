\name{clust_sum}
\alias{clust_sum}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  An intermediate function that summarizes the cl_sample to P.
}
\description{
  P is a list of matrices, with P[[i]] a matrix contain all information of ith cluster assignment. P[[i]] is K by (n+1) matrix. Each row for one cluster, last column as number of targets in the cluster.
}
\usage{
clust_sum(cl_sample, dat, iter, a, b)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cl_sample}{
  A list that contains all the clustering assignment information
}
  \item{dat}{
  Original interaction data x0
}
  \item{iter}{
  Sample size of DPMM clustering samples
}
  \item{a}{
  Prior hyperparameters alpha0
}
  \item{b}{
  Prior hyperparameters beta0
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
#suppose you first have the dat and cl_sample as in dpmm_beta
P = clust_sum(cl_sample, dat, iter = length(cl_sample$KK), a, b)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

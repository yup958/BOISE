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
Update_beta(cl, dat, a, b, alpha = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cl}{
  Old clustering assignment list of {K, N, C}
}
  \item{dat}{
  Original interaction data x0
}
  \item{a}{
  Hyperparameters alpha0
}
  \item{b}{
  Hyperparameters beta0
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
data(pkis1)
## Use 2sd criteria to create binary matrix
dat <- t(apply(pkis1, 1, function(x){
  thres = mean(x) + 2 * sd(x)
  return(as.numeric(x>thres))
}))
rm(pkis1)
a = rep(mean(dat,na.rm = T),ncol(dat))
b = 1-a
alpha = 10
cl = Initial_beta(dat,a,b,alpha)
cl = Update_beta(cl, dat, a, b, alpha)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

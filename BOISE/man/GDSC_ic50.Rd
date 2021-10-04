\name{GDSC_ic50}
\alias{GDSC_ic50}
\docType{data}
\title{
  A complete version of Genomics of Drug Sensitivity in Cancer (GDSC) data. 
}
\description{
281 cancer cell lines, 207 anti-cancer drugs. Complete but not binarized data. 
}
\usage{data("GDSC_ic50")}
\format{
  The format is:
 num [1:281, 1:207] 0.4786 -1.49256 -1.66977 0.00562 -1.18833 ...
 - attr(*, "dimnames")=List of 2
  ..$ : chr [1:281] "MC-CAR" "ES3" "ES5" "ES7" ...
  ..$ : chr [1:207] "BMS-536924" "GSK269962A" "Doxorubicin" "Etoposide" ...
}
\details{
A data matrix consists of interaction between cancer cell lines and anti-cancer drugs
}
\source{
Yang, W., Soares, J., Greninger, P. et al. (2013) Genomics of Drug Sensitivity in Cancer (GDSC): a resource for therapeutic biomarker discovery in cancer cells. Nucleic Acids Research, vol. 41, pp. 955-961. 
}
\references{
Yu, P., Ericksen, S.S., Gitter, A. and Newton, M.A. (2020) \emph{Bayes Optimal Informer Sets for Drug Discovery}. University of Wisconsin, Technical Report, 2020.

}
\examples{
data(GDSC_ic50)
## maybe str(GDSC_ic50) ; plot(GDSC_ic50) ...
}
\keyword{datasets}

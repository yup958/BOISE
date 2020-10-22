\name{GDSC_ic50}
\alias{GDSC_ic50}
\docType{data}
\title{
  GDSC cleaned data.
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
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(GDSC_ic50)
## maybe str(GDSC_ic50) ; plot(GDSC_ic50) ...
}
\keyword{datasets}

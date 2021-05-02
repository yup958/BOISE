\name{pkis1}
\alias{pkis1}
\docType{data}
\title{
  PKIS1 data: GSK Published Kinase Inhibitor Set, a chemogenomic data of inhibition rates between human kinases and inhibitors.
}
\description{
 224 targets and 366 compounds. Complete non-binarized data.
}
\usage{data("pkis1")}
\format{
  The format is:
 num [1:224, 1:366] 0.1397 0.0568 0.1041 0.1502 0.0924 ...
 - attr(*, "dimnames")=List of 2
  ..$ : chr [1:224] "ABL1" "ABL1_E255K" "ABL1_H396P" "ABL1_M351T" ...
  ..$ : chr [1:366] "6539592" "6539600" "6539599" "6539597" ...
}
\details{
A data matrix consists of inhibition rates between human kinases and inhibitors.
}
\source{
H Drewry, D., M Willson, T., & J Zuercher, W. (2014). Seeding collaborations to advance kinase science with the GSK Published Kinase Inhibitor Set (PKIS). Current topics in medicinal chemistry, 14(3), 340-342.
}
\references{

 Zhang H, Ericksen SS, Lee C-p, Ananiev GE, Wlodarchak N, Yu P, et al. (2019) Predicting kinase inhibitors using bioactivity matrix derived informer sets. PLoS Comput Biol 15(8): e1006813. https://doi.org/10.1371/journal.pcbi.1006813
}
\examples{
data(pkis1)
## maybe str(pkis1) ; plot(pkis1) ...
}
\keyword{datasets}

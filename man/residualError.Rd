\name{residualError}
\alias{residualError}
\title{
Compute the residual error for a dataset
}
\description{
This function takes in a CNOlist and computes the residual error, which is the minimum error between the scaled continuous data and a binary boolean approximation of this data.
}
\usage{
residualError(CNOlist)
}
\arguments{
  \item{CNOlist}{
a CNOlist
}
}
\details{
Be aware that it is expected that \code{$valueSignals[[1]]} holds t0 (all signals=0) and \code{$valueSignals[[2]]} holds t1, \code{$valueSignals[[3]]} holds t2. If you give a CNOlist with more than 3 elements in valueSignals you will get a warning message but the function should still work based on the first 2 time points. If you give a CNOlist with only 2 elements in valueSignals (i.e. you don't have a time 2), the function will fill in the residual error t1 and leave t2 and t1andt2 = NA
}
\value{
a vector with named entries t1, t2 and t1andt2 that hold the residual error for when only t1 is considered, only t2 is considered, or both are considered
}
\author{
C. Terfve
}
\seealso{
makeCNOlist, normaliseCNOlist, GetFit
}
\examples{
data(CNOlistToy,package="CellNOptR")
resECNOlistToy<-residualError(CNOlistToy)
}

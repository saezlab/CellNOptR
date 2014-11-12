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
Be aware that it is expected that \code{$valueSignals[[1]]} holds t0 (all
signals=0) and \code{$valueSignals[[2]]} holds t1, \code{$valueSignals[[3]]}
holds t2 and so on. The output is a list of residual errors at each time greater
than 2. In addition, the total error is stored. For back compatibility, an
additional field called t1andt2 is stored (NA is only t1 provided). 

}
\value{
a vector with named entries t1, t2, ...tn, t1andt2 and total. If only t1 is
provided, t1andt2 is set to NA. The field t1andt2 may be removed in the future.
Use the field called total instead.
}
\author{
C. Terfve, T. Cokelaer
}
\seealso{
\link{makeCNOlist}, \link{normaliseCNOlist}, \link{getFit}
}
\examples{
data(CNOlistToy,package="CellNOptR")
resECNOlistToy<-residualError(CNOlistToy)
}

\name{crossInhibitedData}
\alias{crossInhibitedData}
\title{
If an inhibitor is also a measured species, replace the data with NA (when inhibited)
}
\description{
If an inhibitor is also a measured species, replace the data with NA (when inhibited)
}
\usage{
crossInhibitedData(object)
}

\arguments{
  \item{object}{
the CNOlist that contains the data
}
}



\value{
  the new cnolist
}

\author{
 T. Cokelaer
}



\examples{
data(ToyModel,package="CellNOptR")
data(CNOlistToy,package="CellNOptR")
cnolist = crossInhibitedData(CNOlist(CNOlistToy))
}

\name{crossInhibitedData}
\alias{crossInhibitedData}
\title{
If an inhibitor is also a measured species, replace the data with NA (when inhibited)
}
\description{
If an inhibitor is also a measured species, replace the data with NA (when inhibited)
}
\usage{
crossInhibitedData(cnolist)
}

\arguments{
  \item{data}{
the CNOlist that contains the data
}
}



\value{
  the new cnolist
}

\author{
 T. Cokelaer
}


\seealso{

}
\examples{
data(ToyModel,package="CellNOptR")
data(CNOlistToy,package="CellNOptR")
cnolist = crossInhibitedData(CNOlistToy)
}

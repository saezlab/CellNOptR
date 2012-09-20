\name{prep4sim}
\alias{prep4sim}
\alias{prep4Sim}
\title{
Prepare a model for simulation
}
\description{
Adds to the model some fields that are used by the simulation engine
}
\usage{
prep4sim(model)
}
\arguments{
  \item{model}{
a model list, as output by readSIF, normally pre-processed but that is not a requirement of this function
}
}
\details{
This adds fields that are necessary for the simulation engine in a version that is extensible for constrained Fuzzy logic extension of the methods applied here (in development).
}
\value{
this function returns a list with fields:
\item{finalCube}{stores, for each reac(row) the location of its inputs (col)}
\item{ixNeg}{stores, for each reac(row) and each input (col) whether it is a negative input}
\item{ignoreCube}{logical matrix of the same size as the 2 above, that tells whether the particular cell is filled or not}
\item{maxIx}{row vector that stores, for each reac, the location of its output}
\item{modelname}{stores the name of the model from which these fields were derived}
\item{maxInput}{stores the max number of inputs observed in the model for a single reaction}
}
\author{
C. Terfve, T. Cokelaer
}
\seealso{
\link{simulatorT1}
}
\note{For internal usage since version 1.3.28}
\examples{
data(ToyModel,package="CellNOptR")
ToyFields4Sim<-prep4sim(ToyModel)
}

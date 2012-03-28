\name{simulateT1}
\alias{simulateT1}
\title{
Cut and simulation of a boolean model at t1
}
\description{
This function cuts a model according to a bitstring optimised at T1, and
simulates the model accordingly.
}
\usage{
simulateT1(CNOlist, Model, bStringT1, SimList, indexList)
}
\arguments{
  \item{CNOlist}{
a CNOlist object
}
  \item{Model}{
a full model
}
  \item{bStringT1}{
a bitstring to cut the model, as output by \code{gabinaryT1} (i.e. a vector of
1s and 0s, of length equal to the number of reactions in the model)
}
  \item{SimList}{
a list of additional fields for simulation as created by \code{prep4Sim},
corresponding to the full model
}
  \item{indexList}{
a list of indexes as created by indexFinder
}
}
\details{
This function is a wrapper for simulatorT1, that cuts the model before simulating it
}
\value{
a matrix of simulated values, including all species in the model, i.e. to be
used as input of gabinaryT2 (not implemented here) but not to be used directly
in \code{plotOptimResults.}
}
\author{
C.Terfve
}
\seealso{
cutAndPlotOptimResultsT1, simulatorT1
}
\examples{
# This will compute the output of a random model obtained by randomly selecting
# which gates of the initial models are included.

data(CNOlistToy,package="CellNOptR")
data(ToyModel,package="CellNOptR")

indicesToy<-indexFinder(CNOlistToy,ToyModel,verbose=FALSE)
ToyFields4Sim<-prep4Sim(ToyModel)

SimRes<-simulateT1(
	CNOlist=CNOlistToy, 
	Model=ToyModel, 
	bStringT1=round(runif(length(ToyModel$reacID))), 
	SimList=ToyFields4Sim, 
	indexList=indicesToy)
}

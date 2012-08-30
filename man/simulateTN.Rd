\name{simulateTN}
\alias{simulateTN}
\title{
Cut and simulation of a boolean model at t1
}
\description{
This function cuts a model according to a bitstring optimised at T1, and
simulates the model accordingly.
}
\usage{
simulateTN(CNOlist, model, bStrings)
}
\arguments{
  \item{CNOlist}{
a CNOlist object
}
  \item{model}{
a full model
}
  \item{bStrings}{
a bitstring to cut the model, as output by \code{gaBinaryT1} (i.e. a vector of
1s and 0s, of length equal to the number of reactions in the model)
}
}
\details{
This function is a wrapper around the family of functions called simulatorT1 ,
T2 and TN.
}
\value{
a matrix of simulated values, including all species in the model, i.e. to be
used as input of gaBinaryTN.
}
\author{
T.Cokelaer, S.Schrier based on simulatorT1 (C.Terfve)
}
\seealso{
\link{cutAndPlotResultsT1}, \link{simulatorT1}
}
\examples{
# This will compute the output of a random model obtained by randomly selecting
# which gates of the initial models are included.

data(CNOlistToy,package="CellNOptR")
data(ToyModel,package="CellNOptR")

simRes<-simulateTN(
	CNOlist=CNOlistToy, 
	model=ToyModel, 
	bStrings=list(round(runif(length(ToyModel$reacID)))))
}

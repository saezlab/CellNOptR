\name{simulatorT2}
\alias{simulatorT2}
\title{
Deprecated since 1.3.28. Use \link{simulateTN} function instead.
Simulation of a boolean model for time 2
}
\description{
This function simulates a boolean model at time 2 where time 2 is assume to be a pseudo-steady states at a time scale slower than the pseudo-steady state evaluated at time 1
}
\usage{
simulatorT2(simResultsT1, CNOlist, model, simList, indexList, timeIndex=3)
}
\arguments{
  \item{simResultsT1}{
a matrix that is the output of simulatorT1 (i.e. one row per condition and one column per species IN THE MODEL)
}
  \item{CNOlist}{
a CNOlist
}
  \item{model}{
a Model that only contains the reactions to be evaluated, and the additional field \code{model$times} that should have been created inside the gabinaryT2 optimisation engine.}
  \item{simList}{
a simList as created by prep4sim, that has also already been cut to contain only the reactions to be evaluated
}
  \item{indexList}{
an indexList as created by indexFinder
}
\item{timeIndex}{
argument requiter only for steady states T3 and above to specify the Time
indices.
}
}
\details{
This is the simulator for time T2, it is very similar to simulatorT1 but here we assume that we start from the simulated results at t1 (i.e. we start from a pseudo-steady state) then it does a first iteration, and whatever branch is set to be active at t2 has an effect that cannot be changed after the first iteration, i.e. the output node of a t2 iteration is fixed. We assume here that the model has already been cut, and that the cutting is based on keeping  all the edges that are set to either 1 or 2, and there is an additional field \code{$times} in the model that keeps the info of the t1/t2 (it is a vector of 1s and 2s of length=number of reaches present). Structurally the function is almost identical to simulatorT1 but it does a first iteration where all the gates that lead to a node that also receives a T2 gates are set to the same value as the t2 gate in the ANDs calculation.In the main loop, the nodes that are targets of t2 interactions are constantly reset to their value at the first iteration, at the end of each iteration (similarly to what is done with stimulated and inhibited species)

The \code{model$times} field is a vector of 1s and 2s that tells the simulator which interactions are expected to be active at t1 and which are at t2.
}
\value{
This function outputs a single matrix of format similar to valueSignals in the CNOlist but that contains an output for each species in the model.  This matrix is the simulated equivalent of valueSignals at time 2 if you consider only the columns given by indexSignals
}
\author{
C. Terfve
}
\seealso{
\link{simulateTN}, \link{cutAndPlotResultsT1}, \link{simulatorT1}
}
\examples{

# This computes the output of the full model, which is normally not done on a
# stand alone basis, but if you have a model and would like to visualise its
# output compared to your data, then this is what you should do

data(CNOlistToy2,package="CellNOptR")
data(ToyModel2,package="CellNOptR")

#Sim<-simulatorTN(
#   CNOlist=CNOlistToy2,
#   model=ToyModel2)

#Sim2<-simulatorTN(
#   simResultsT1=Sim, 
#  CNOlist=CNOlistToy2, 
#   model=ToyModel2, 
#   simList=ToyFields4Sim, 
#   indexList=indicesToy, timeIndex=3)


}


\name{getFit}
\alias{getFit}
\title{
Compute the score of a model
}
\description{
This function computes the value of the objective function for a model and an associated data set, as a sum of a term that computes the fit of model to data, a term that penalises the NA values produced by the model, and a term that penalises increasing size of the model.
}
\usage{
getFit(SimResults, CNOlist, Model, indexList, timePoint = c("t1", "t2"), sizeFac = 1e-04, NAFac = 1)
}
\arguments{
  \item{SimResults}{
matrix of simulated results (the full one as output by the simulator)
}
  \item{CNOlist}{
a CNOlist to compare the simulated results with
}
  \item{Model}{
a model that has already been cut to contain only the reactions in the optimal bitstring
}
  \item{indexList}{
list of indexes as produced by indexFinder
}
  \item{timePoint}{
"t1" or "t2" tells which time point we are looking at. If timePoint=t1 then we will compare the SimResults to the results stored in CNOlist$valueSignals[[2]].  If timePoint=t1 then we will compare the SimResults to the results stored in CNOlist$valueSignals[[2]]
}
  \item{sizeFac}{
weights the penalty for the size of the model, default=0.0001
}
  \item{NAFac}{
weights the penalty for the number of NAs
}
}
\details{
BE AWARE: contrary to what is done in the Matlab version of CellNOpt, here the simulation results are computed beforehand and the Model that is input into this function is a model that has already been cut i.e. that only contains the reactions present in the optimised model (i.e.should be the same model as the one that you input into the simulator). Also, the SimResults matrix is the full one as output by the simulator, i.e. it contains results for all species in the model, not only the signals
}
\value{
This function returns a single number, the value of the objective function.
}
\references{
J. Saez-Rodriguez, L. G. Alexopoulos, J. Epperlein, R. Samaga, D. A. Lauffenburger, S. Klamt and P. K. Sorger. Discrete logic modeling as a means to link protein signaling networks with functional analysis of
mammalian signal transduction, Molecular Systems Biology, 5:331, 2009.
}
\author{
C. Terfve
}
\seealso{
gabinaryT1, indexFinder, simulatorT1
}
\examples{
#Here we will evaluate the fit of the full initial model, 
#without pre-processing or any optimisation

data(CNOlistToy,package="CellNOptR")
data(ToyModel,package="CellNOptR")

indicesToy<-indexFinder(CNOlistToy,ToyModel,verbose=FALSE)
ToyFields4Sim<-prep4Sim(ToyModel)
SimResults<-simulatorT1(
	CNOlist=CNOlistToy,
	Model=ToyModel,
	SimList=ToyFields4Sim,
	indexList=indicesToy)
Score<-getFit(
	SimResults=SimResults,
	CNOlist=CNOlistToy,
	Model=ToyModel,
	indexList=indicesToy,
	timePoint="t1")
}

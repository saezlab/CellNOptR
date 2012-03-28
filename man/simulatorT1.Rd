\name{simulatorT1}
\alias{simulatorT1}
\title{
Simulation of a boolean model
}
\description{
This is the simulator, inspired from BoolSimEngMKM in the Matlab CellNOpt, to be used on one time point simulations
}
\usage{
simulatorT1(CNOlist, Model, SimList, indexList)
}
\arguments{
  \item{CNOlist}{
a CNOlist
}
  \item{Model}{
a Model that only contains the reactions to be evaluated
}
  \item{SimList}{
a SimList as created by prep4Sim, that has also already been cut to contain only the reactions to be evaluated
}
  \item{indexList}{
an indexList as created by indexFinder
}
}
\details{
Differences from the BoolSimEngMKM simulator include: 
the valueInhibitors has not been previously flipped;
the function outputs the values across all conditions for all species in the model, instead of only for the signal species. This is because then the  output of this function can be used as initial values for the version of the simulator that works on time point 2 (not implemented in this version).

If you would like to compute the output of a model that contains some of the gates in the model but not all, we suggest that you use the function \code{SimulateT1} and specify in the \code{bStringT1} argument which gates you want to be included.  Indeed, SimulateT1 is a wrapper around simulatorT1 that takes care of cutting the model for you before simulating it.
}
\value{
This function outputs a single matrix of format similar to valueSignals in the CNOlist but that contains an output for each species in the model.  This matrix is the simulated equivalent of valueSignals at time 1, if you consider only the columns given by indexSignals.
}
\references{
\enumerate{
\item{J. Saez-Rodriguez, L. G. Alexopoulos, J. Epperlein, R. Samaga, D. A.
Lauffenburger, S. Klamt and P. K. Sorger. Discrete logic modeling as a means to
link protein signaling networks with functional analysis of mammalian signal
transduction, Molecular Systems Biology, 5:331, 2009.}
\item{M. K. Morris, J. Saez-Rodriguez, D. Clarke, P. K. Sorger, D. A.
Lauffenburger. Training Signaling Pathway Maps to Biochemical Data with
Constrained Fuzzy Logic: Quantitative Analysis of Liver Cell Responses to
Inflammatory Stimuli, PLoS Comp. Biol., 7(3): e1001099, 2011.}
}
}
\author{
C. Terfve
}
\seealso{
SimulateT1, cutAndPlotResultsT1
}
\examples{
#This computes the output of the full model, which is normally not done on a stand alone basis, but if you have a model and would like to visualise its output compared to your data, then this is what you should do

data(CNOlistToy,package="CellNOptR")
data(ToyModel,package="CellNOptR")

indicesToy<-indexFinder(CNOlistToy,ToyModel,verbose=TRUE)
ToyFields4Sim<-prep4Sim(ToyModel)

Sim<-simulatorT1(
	CNOlist=CNOlistToy,
	Model=ToyModel,
	SimList=ToyFields4Sim,
	indexList=indicesToy)
}

\name{cSimulator}
\alias{cSimulator}
\title{
C-implementation of simulatorT1.
}
\description{
This is the simulator, inspired from BoolSimEngMKM in the Matlab CellNOpt, to be
used on one time point simulations. Use the R interface provided by
\link{simulatorT1}.
}
\usage{
cSimulator(CNOlist, model, simList, indexList)
}
\arguments{
  \item{CNOlist}{
a CNOlist
}
  \item{model}{
a model that only contains the reactions to be evaluated
}
  \item{simList}{
a simList as created by prep4sim, that has also already been cut to contain only the reactions to be evaluated
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
A. MacNamara based on former simulatorT1 version from C.Terfve.
}
\seealso{
\link{simulatorT1}.
}

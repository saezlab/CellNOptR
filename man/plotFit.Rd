\name{plotFit}
\alias{plotFit}
\title{
Plot the evolution of an optimisation
}
\description{
This function takes in the results of an optimisation by gaBinaryT1 and plots the evolution of best fit and average fit against generations.
}
\usage{
plotFit(OptRes)
}
\arguments{
  \item{OptRes}{
an object created by the optimisation engine (gabinaryT1)
}
}
\value{
This function doesn't return anything, it just produces a plot in your graphics window.
}
\author{
C. Terfve
}
\seealso{
gabinaryT1
}
\examples{
data(CNOlistToy,package="CellNOptR")
data(ToyModel,package="CellNOptR")

#process the model

indicesToy<-indexFinder(CNOlistToy,ToyModel,verbose=FALSE)
ToyComp<-compressModel(ToyModel,indicesToy)
indicesToyComp<-indexFinder(CNOlistToy,ToyComp)
ToyCompExp<-expandGates(ToyComp)

#optimise

ToyFields4Sim<-prep4Sim(ToyCompExp)
initBstring<-rep(1,length(ToyCompExp$reacID))
ToyT1opt<-gaBinaryT1(
	CNOlist=CNOlistToy,
	Model=ToyCompExp,
	SimList=ToyFields4Sim,
	indexList=indicesToyComp,
	initBstring=initBstring,
	maxGens=2,
	PopSize=5,
	verbose=FALSE)

plotFit(OptRes=ToyT1opt)
}

\name{cutAndPlotResultsT2}
\alias{cutAndPlotResultsT2}
\title{
Plot the results of an optimisation at t2
}
\description{
This function takes a model and an optimised bitstring, it cuts the model according to the bitstring and plots the results of the simulation along with the experimental data. This function is designed to work on results of a 2 step optimisation.
}
\usage{
cutAndPlotResultsT2(Model, bStringT1, bStringT2, SimList, CNOlist, indexList, plotPDF = FALSE, tag=NULL)
}
\arguments{
  \item{Model}{
a model (the full one that was used for optimisation)
}
  \item{bStringT1}{
a bitstring for T1 as output by gabinaryT1 
}
  \item{bStringT2}{
a bitstring for T2 as output by gabinaryT2
}
  \item{SimList}{
a simlist corresponding to the model, as output by prep4Sim
}
  \item{CNOlist}{
a CNOlist, corresponding to the optimisation one
}
  \item{indexList}{
an indexList, produced by indexFinder ran on the model and the CNOlist above
}
  \item{plotPDF}{
TRUE or FALSE, tells whether you want a pdf to be produced or not
}
  \item{tag}{
NULL or string; tells whether you want to prefix filenames with a tag (replaces
the default behaviour).
}
}

\value{
This function doesn't return anything, it only plots the graph in your graphic window, and in a pdf file if asked
}
\author{
C. Terfve
}
\seealso{
gabinaryT1, prep4Sim, cutAndPlotResultsT1
}
\examples{
#load data

data(CNOlistToy2,package="CellNOptR")
data(ToyModel2,package="CellNOptR")

#pre-process model

checkSignals(CNOlistToy2,ToyModel2)
indexesToy2<-indexFinder(CNOlistToy2,ToyModel2,verbose=FALSE)
ToyNCNOindexes2<-findNONC(ToyModel2,indexesToy2,verbose=FALSE)
ToyNCNOcut2<-cutNONC(ToyModel2,ToyNCNOindexes2)
indexesToyNCNOcut2<-indexFinder(CNOlistToy2,ToyNCNOcut2)
ToyNCNOcutComp2<-compressModel(ToyNCNOcut2,indexesToyNCNOcut2)
indexesToyNCNOcutComp2<-indexFinder(CNOlistToy2,ToyNCNOcutComp2)
ToyNCNOcutCompExp2<-expandGates(ToyNCNOcutComp2)

#optimise t1
ToyFields4Sim2<-prep4Sim(ToyNCNOcutCompExp2)
initBstring2<-rep(1,length(ToyNCNOcutCompExp2$reacID))
ToyT1opt2<-gaBinaryT1(
	CNOlist=CNOlistToy2,
	Model=ToyNCNOcutCompExp2,
	SimList=ToyFields4Sim2,
	indexList=indexesToyNCNOcutComp2,
	initBstring=initBstring2,
	maxGens=2, 
	PopSize = 5,
	verbose=FALSE)

#Optimise T2
SimToyT12<-simulateT1(
	CNOlist=CNOlistToy2,
	Model=ToyNCNOcutCompExp2,
	bStringT1=ToyT1opt2$bString,
	SimList=ToyFields4Sim2,
	indexList=indexesToyNCNOcutComp2)
ToyT1opt2T2<-gaBinaryT2(
	CNOlist=CNOlistToy2,
	Model=ToyNCNOcutCompExp2,
	SimList=ToyFields4Sim2,
	indexList=indexesToyNCNOcutComp2,
	bStringT1=ToyT1opt2$bString,
	SimResT1=SimToyT12,
	maxGens=2, 
	PopSize = 5,
	verbose=FALSE)
cutAndPlotResultsT2(
	Model=ToyNCNOcutCompExp2,
	bStringT1=ToyT1opt2$bString,
	bStringT2=ToyT1opt2T2$bString,
	SimList=ToyFields4Sim2,
	CNOlist=CNOlistToy2,
	indexList=indexesToyNCNOcutComp2,
	plotPDF=FALSE)
}

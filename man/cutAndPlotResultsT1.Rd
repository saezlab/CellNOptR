\name{cutAndPlotResultsT1}
\alias{cutAndPlotResultsT1}
\title{
Plot the results of an optimisation at t1
}
\description{
This function takes a model and an optimised bitstring, it cuts the model according to the bitstring and plots the results of the simulation along with the experimental data.
}
\usage{
cutAndPlotResultsT1(Model, bString, SimList, CNOlist, indexList, plotPDF =
FALSE, tag = NULL, show=TRUE)
}
\arguments{
  \item{Model}{
a model (the full one that was used for optimisation)
}
  \item{bString}{
a bitstring for T1  as output by gabinaryT1 (i.e. a vector of 1s and 0s)
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
TRUE or FALSE; tells whether you want a pdf to be produced or not
}
  \item{tag}{
NULL or string; tells whether you want to prefix filenames with a tag (replaces
the default behaviour).
}
  \item{show}{
TRUE of FALSE; FALSE tells not to show the plot, which can be useful when scripting.
}
}

\value{
This function doesn't return anything, it only plots the graph in your graphic window, and in a pdf file if asked
}
\author{
C.Terfve
}
\seealso{
gabinaryT1, prep4Sim
}
\examples{
#load data

data(CNOlistToy,package="CellNOptR")
data(ToyModel,package="CellNOptR")

#pre-process model

checkSignals(CNOlistToy,ToyModel)
indicesToy<-indexFinder(CNOlistToy,ToyModel,verbose=FALSE)
ToyNCNOindices<-findNONC(ToyModel,indicesToy,verbose=FALSE)
ToyNCNOcut<-cutNONC(ToyModel,ToyNCNOindices)
indicesToyNCNOcut<-indexFinder(CNOlistToy,ToyNCNOcut)
ToyNCNOcutComp<-compressModel(ToyNCNOcut,indicesToyNCNOcut)
indicesToyNCNOcutComp<-indexFinder(CNOlistToy,ToyNCNOcutComp)
ToyNCNOcutCompExp<-expandGates(ToyNCNOcutComp)

#optimise

ToyFields4Sim<-prep4Sim(ToyNCNOcutCompExp)
initBstring<-rep(1,length(ToyNCNOcutCompExp$reacID))
ToyT1opt<-gaBinaryT1(
	CNOlist=CNOlistToy,
	Model=ToyNCNOcutCompExp,
	SimList=ToyFields4Sim,
	indexList=indicesToyNCNOcutComp,
	initBstring=initBstring, 
	maxGens=2, 
	PopSize = 5,
	verbose=FALSE)
cutAndPlotResultsT1(
	Model=ToyNCNOcutCompExp,
	bString=ToyT1opt$bString,
	SimList=ToyFields4Sim,
	CNOlist=CNOlistToy,
	indexList=indicesToyNCNOcutComp,
	plotPDF=FALSE)
}

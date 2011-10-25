\name{writeNetwork}
\alias{writeNetwork}
\title{
Write a previous knowledge network model to a sif file (with attribute files), as well as a dot file
}
\description{
This function writes the original previous knowledge network (the model that you loaded in the beginning of your analysis)  in a sif file, with a nodes attribute file that specifies if each node was stimulated/inhibited/signal/compressed/non-controllable-non-observable and an edge attribute file that specifies if the edge was absent in the optimal model (0) present in the optimal model at t1 (1) or present in the optimal model at t2 (2).

This function also writes a Graphviz dot file that contains the same information (see \code{writeDot} for more information about the dot file conventions).
}
\usage{
writeNetwork(ModelOriginal, ModelComprExpanded, optimResT1, optimResT2, CNOlist)
}
\arguments{
  \item{ModelOriginal}{
The PKN model
}
  \item{ModelComprExpanded}{
The scaffold model (i.e. compressed and expanded)
}
  \item{optimResT1}{
The results of the optimisation process at t1
}
  \item{optimResT2}{
The results of the optimisation process at t2 (set this to NA, the t2 optimisation is not implemented in this version).
}
  \item{CNOlist}{
The CNOlist on which the optimisation is based
}
}
\details{
The weights of the edges are computed as the mean across models within the relative tolerance limits, as output in the results from the optimisation \code{$StringsTol}. Strings that are in \code{$StringsTol} are the ones that are within the relative tolerance limits around the best solution in the population across all generations of the optimisation. 

!If there is no time 2, then the argument optimResT2 should be = NA

This function maps back the edges weights from the optimised (expanded and compressed) model to the original model. The mapping back only works if the path has length 2 at most (i.e. you have node1-comp1-comp2-node2, where comp refer to nodes that have been compressed).
}
\value{
This function does not have any output, it just writes a sif file, an edge attribute file, and a node attribute file
}
\author{
C. Terfve
}
\note{
The mapback of this function is still an open question, even in the Matlab version.  Future developments will include more robust versions of the mapping back algorithm, probably as a separate mapback function.
}
\seealso{
writeScaffold, writeDot
}
\examples{
tmpdir<-tempdir()
setwd(tmpdir)

#load data

data(CNOlistToy,package="CellNOptR")
data(ToyModel,package="CellNOptR")

#pre-process model

indicesToy<-indexFinder(CNOlistToy,ToyModel,verbose=TRUE)
ToyNCNOindices<-findNONC(ToyModel,indicesToy,verbose=TRUE)
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
	verbose=TRUE,
	maxGens=2,
	PopSize=5)

#write network

writeNetwork(
	ModelOriginal=ToyModel,
	ModelComprExpanded=ToyNCNOcutCompExp,
	optimResT1=ToyT1opt,
	optimResT2=NA,
	CNOlist=CNOlistToy)
}

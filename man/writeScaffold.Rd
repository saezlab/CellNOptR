\name{writeScaffold}
\alias{writeScaffold}
\title{
Writes the scaffold network to a sif file (with attributes) and to a dot file
}
\description{
This function writes a cytoscape SIF file for the scaffold network, with an associated edge attribute file that holds whether the edge is present at t1,t2 or not present at all and another associated edge attribute file that holds the weights of the edges. This function also writes a dot file that contains the same information (see \code{writeDot} for more information about the dot file conventions).
}
\usage{
writeScaffold(modelComprExpanded, optimResT1, optimResT2, modelOriginal,
CNOlist, tag=NULL)
}
\arguments{
  \item{modelComprExpanded}{
The scaffold model (i.e. compressed and expanded)
}
  \item{optimResT1}{
The results of the optimisation process at t1
}
  \item{optimResT2}{
The results of the optimisation process at t2 (set this to NA if you have performed a one time point optimisation).
}
  \item{modelOriginal}{
The PKN model
}
  \item{CNOlist}{
The CNOlist on which the optimisation is based 
}
  \item{tag}{
NULL or string; tells whether you want to prefix filenames with a tag (replaces the default behaviour).
}
}
\details{
By scaffold network we mean the network that is used as a basis for optimisation (i.e. a compressed and expanded network), therefore no map back of the weights is necessary here.

The weights of the edges are computed as the mean across models within the relative tolerance limits, as output in the results from the optimisation \code{$stringsTol}. Strings that are in \code{$stringsTol} are the ones that are within the relative tolerance limits around the best solution in the population across all generations of the optimisation. 

!If there is no time 2, then the argument optimResT2 should be = NA.

}
\value{
This function does not return anything, it writes a sif file and 2 edge attributes files, and a dot file, in your working directory.
}
\author{
C.Terfve
}
\seealso{
writeNetwork, writeDot
}
\examples{
tmpdir<-tempdir()
setwd(tmpdir)

#load the data

data(CNOlistToy,package="CellNOptR")
data(ToyModel,package="CellNOptR")

#pre-process the model (partial)

indicesToy<-indexFinder(CNOlistToy,ToyModel,verbose=TRUE)
ToyNCNOcutComp<-compressModel(ToyModel,indicesToy)
indicesToyNCNOcutComp<-indexFinder(CNOlistToy,ToyNCNOcutComp)
ToyNCNOcutCompExp<-expandGates(ToyNCNOcutComp)

#optimise

ToyFields4Sim<-prep4sim(ToyNCNOcutCompExp)
initBstring<-rep(1,length(ToyNCNOcutCompExp$reacID))
ToyT1opt<-gaBinaryT1(
	CNOlist=CNOlistToy,
	model=ToyNCNOcutCompExp,
	initBstring=initBstring, 
	maxGens=3,
	popSize=5,
	verbose=TRUE)

#write the network

writeScaffold(
	modelOriginal=ToyModel,
	modelComprExpanded=ToyNCNOcutCompExp,
	optimResT1=ToyT1opt,
	optimResT2=NA,
	CNOlist=CNOlistToy)
}

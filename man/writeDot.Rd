\name{writeDot}
\alias{writeDot}
\title{
Write a model, and attached features, to a dot file 
}
\description{
This function writes a model to a Graphviz dot file with encoded features such as edge weight and nodes status (see details).
}
\usage{
writeDot(dotNodes,dotMatrix,Model,filename)
}
\arguments{
  \item{dotNodes}{
internal variables created by writeNetwork or writeScaffold; dotNodes is a matrix with 2 columns: the first has the node names,and the second the attributes (signal, stimulated, inhibited, compressed, nano). A node can appear twice in this matrix if it belongs to more of one of the above categories; a node could also not appear here if it is is none of these categories
}
  \item{dotMatrix}{
internal variables created by writeNetwork or writeScaffold; dotMatrix is a matrix with 4 or 5 columns, and a row for each reaction:the first column holds the name of the input node, the second column holds the sign of the reaction (-1 if negative, 1 if positive), the third column holds the name of the output node, the fourth column holds the time stamp (0,1,2), an optional 5th column holds the weights of the edges
}
  \item{Model}{
A model to be plotted, if used inside writeNetwork then this should be the previous knowledge network (ModelOriginal), if inside writeScaffold then this should be the scaffold (ModelComprExpanded)
}
  \item{filename}{
a name for the file
}
}
\details{
This function is not to be used on its own, it should be used internally to writeNetwork or writeScaffold. For the colouring of the nodes, nodes that are both stimulated and inhibited or any other combination, only one colour per category is used, and the following order of priority for the colours is used: signals prime over inhibited nodes which primes over stimulated nodes which primes over non-controllable/non-observable nodes, which primes over compressed.  Nodes that are neither of those have a black contour, stimulated nodes are green, inhibited are red, measure are blue, compressed and non-controllable/non-observable nodes are black and dashed.  Edges are coloured according to time stamp in the optimal model (green=t1, blue=t1 and/or t2, grey=neither); on the scaffold, the strokes of the edges reflects the weights in the models within reltol (i.e. for each edge, the weight is the frequency with which it appeared among the models within the relative tolerance boundaries around the best solution).
}
\value{
This function does not have any output, it just writes a dot file in your working directory. 
}
\references{
Emden R. Gansner ,  Stephen C. North. An Open Graph Visualization System and Its Applications to Software Engineering. Software - Practice and Experience (1999)
}
\author{
C. Terfve
}
\seealso{
writeNetwork, writeScaffold
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
	maxGens=2,
	PopSize=5, 
	verbose=TRUE)

#write network

writeNetwork(
	ModelOriginal=ToyModel,
	ModelComprExpanded=ToyNCNOcutCompExp,
	optimResT1=ToyT1opt,
	optimResT2=NA,
	CNOlist=CNOlistToy)
}

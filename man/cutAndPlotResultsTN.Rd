\name{cutAndPlotResultsTN}
\alias{cutAndPlotResultsTN}
\title{
Plot the results of an optimisation at tN
}
\description{
This function takes a model and an optimised bitstring, it cuts the model according to the bitstring and plots the results of the simulation along with the experimental data. This function is designed to work on results of a 2 step optimisation.
}


\usage{
cutAndPlotResultsTN(CNOlist, model, bStrings, plotPDF = FALSE, tag=NULL, 
    plotParams = list(maxrow = 10))
}
\arguments{
  \item{CNOlist}{
a CNOlist, corresponding to the optimisation one
}
  \item{model}{
a model (the full one that was used for optimisation)
}
  \item{bStrings}{
a list of bitstring at different time points 
}
  \item{plotPDF}{
TRUE or FALSE, tells whether you want a pdf to be produced or not
}
  \item{tag}{
NULL or string; tells whether you want to prefix filenames with a tag (replaces
the default behaviour).
}
  \item{maxrow}{
maximum number of row in the plot. 
}
\item{plotParams}{
a list of option related to the PDF and plotting outputs. (1) maxrow
is the maximum number of row used to plot the results. See  \link{plotOptimResultsPan}
for other fields.}
}
\note{New in version 1.3.28}
\value{
This function doesn't return anything, it only plots the graph in your graphic window, and in a pdf file if asked
}
\author{
T. Cokelaer, A. MacNamara, Sarah Schrier, C. Terfve based on \link{cutAndPlotResultsT1}}

\seealso{
\link{gaBinaryT1}, \link{prep4sim}, \link{cutAndPlotResultsT1}
}
\examples{
#load data

data(CNOlistToy2,package="CellNOptR")
data(ToyModel2,package="CellNOptR")

#pre-process model
model = preprocessing(CNOlistToy2, ToyModel2)

#optimise t1
ToyT1<-gaBinaryT1(
	CNOlist=CNOlistToy2,
	model=model,
	maxGens=2, 
	popSize = 5,
	verbose=FALSE)

#Optimise T2
ToyT2<-gaBinaryTN(
	CNOlist=CNOlistToy2,
	model=model,
	bString=list(ToyT1$bString),
	maxGens=2, 
	popSize = 5,
	verbose=FALSE)

cutAndPlotResultsTN(
	CNOlist=CNOlistToy2,
	model=model,
	bStrings=list(ToyT1$bString, bStringT2=ToyT2$bString),
	plotPDF=FALSE)
}

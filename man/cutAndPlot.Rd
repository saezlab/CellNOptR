\name{cutAndPlot}
\alias{cutAndPlot}
\title{
Interface to cutAndPlotResults functions.
}
\description{
This function takes a model and cnolist as well as a list of optimised bitstring at different time points. It calls the appropriate cutAndPlotResultsTX function.
}
\usage{
cutAndPlot(CNOlist, model,  bStrings, plotPDF=FALSE, tag=NULL,
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
a bitstring for T1  as output by gaBinaryT1 (i.e. a vector of 1s and 0s)
}
  \item{plotPDF}{
TRUE or FALSE; tells whether you want a pdf to be produced or not
}
  \item{tag}{
NULL or string; tells whether you want to prefix filenames with a tag (replaces
the default behaviour).
}
\item{plotParams}{
a list of option related to the PDF and plotting outputs. (1) maxrow
is the maximum number of row used to plot the results. See
\link{plotOptimResultsPan} for other fields.
}

}
\value{
This function returns nothing. It plots a graph in your graphic window and sqve it in a file if asked
}
\author{
T. Cokelaer
}
\seealso{
\link{cutAndPlotResultsT1}
}
\examples{
#load data

data(CNOlistToy,package="CellNOptR")
data(ToyModel,package="CellNOptR")

#pre-process model
model = preprocessing(CNOlistToy, ToyModel)

cutAndPlot(CNOlistToy, model,
    bStrings=list(rep(1,length(model$reacID))),	plotPDF=FALSE)
}

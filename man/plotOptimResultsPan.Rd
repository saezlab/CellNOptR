\name{plotOptimResultsPan}
\alias{plotOptimResultsPan}
\title{
Plots the data and simulated values from any CellNOptR formalism
}
\description{
This function plots the data and simulated values according to each experiment in CNOlist. The data is shown as black triangles and the simulation by a blue dashed line. The combination of cues is given by a panel where black denotes the presence and white the absence of the cue. The goodness-of-fit between model and data is color-coded on a continuous scale from white to red.
}
\usage{
plotOptimResultsPan(simResults, yInterpol=NULL, xCoords=NULL,
CNOlist=CNOlist, formalism=c("ss1","ss2","ssN","dt","ode"), pdf=FALSE,
pdfFileName="", tPt=NULL, 
plotParams = list(margin = 0.1, width=15, height=12, cmap_scale=1, cex=1.6,
ymin=NULL))

}
\arguments{
\item{simResults}{
A list with a field for each time point, each containing a matrix of dimensions (number of conditions) * (number of signals), with the first field being t0.  Typically produced by simulating a model and then extracting the columns that correspond to signals.
}
\item{yInterpol}{
If using CNORdt, these are the interpolated experimental results from getFitTimeScale() that are needed to compare against the Boolean simulation.
}
\item{xCoords}{
These are the x-coordinates obtained from the optimized scaling factor in CNORdt that allow for comparison between time course experimental data and a Boolean model.
}
\item{CNOlist}{
A CNOlist. 
}
\item{formalism}{
An abbreviation of the CellNOptR formalism being used.
}
\item{pdf}{
A Boolean argument denoting whether to print the figure produced by this function to file.
}
\item{pdfFileName}{
If printing to file, the filename to be used.
}
\item{tPt}{
The number of time points in the data.
}


\item{plotParams}{
a list of option related to the PDF and plotting outputs. Currently, the
following attributes are used: (1) margin of the boxes, (2) width and heigth used while creating the
PDF, (3) cmap_scale a value that scales the colors towards small errors (<1) or
large errors (>1); default is 1 (linear colormap) (4) cex is the fontsize used
in the header (5) ymin sets the minimum y axis limit; by default it is the
minimum value found over all data points and therefore can be negative.
}
}
\details{
Depending on the logic formalism, this function is generally called from cutAndPlotResults*(). As shown in the example below however, it can plot the fit of any data and corresponding compatible model. The color denotes the goodness-of-fit, where white shows no difference between simulation and data and red is the maximum error from all conditions and readouts.
}
\value{
This function does not return a value.
}
\references{
J. Saez-Rodriguez, L. G. Alexopoulos, J. Epperlein, R. Samaga, D. A. Lauffenburger, S. Klamt and P. K. Sorger. Discrete logic modeling as a means to link protein signaling networks with functional analysis of
mammalian signal transduction, Molecular Systems Biology, 5:331, 2009.
}
\author{
A. MacNamara
}
\seealso{
cutAndPlotResultsT1
}
\examples{

data(CNOlistToy,package="CellNOptR")
data(ToyModel,package="CellNOptR")
indicesToy <- indexFinder(CNOlistToy, ToyModel, verbose=TRUE)
ToyFields4Sim <- prep4sim(ToyModel)

# simulate model
simRes <- simulatorT1(CNOlist=CNOlistToy, model=ToyModel, simList=ToyFields4Sim, indexList=indicesToy)
simRes = simRes[, indicesToy$signals]

# format data and results

simResults <- list(t0=matrix(data=0, nrow=dim(simRes)[1], ncol=dim(simRes)[2]), t1=simRes)
# plot
plotOptimResultsPan(simResults,
  CNOlist=CNOlistToy,
  formalism="ss1",
  pdf=FALSE,
  tPt=10
)
}

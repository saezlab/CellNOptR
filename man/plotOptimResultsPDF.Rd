\name{plotOptimResultsPDF}
\alias{plotOptimResultsPDF}
\title{
Plot the data and simulated values in a pdf file
}
\description{
This is a wrapper for \code{plotOptimResults}
}
\usage{
plotOptimResultsPDF(SimResults = SimResults, expResults = expResults, times = times, namesCues = namesCues, namesSignals = namesSignals, valueCues = valueCues, filename)
}
\arguments{
  \item{SimResults}{
a list with a field for each time point, each containing a matrix of dimensions number of conditions * number of signals, with the first field being t0.  Typically produced by simulating a model and then extracting the columns that correspond to signals
}
  \item{expResults}{
same as above, but contains the experimental results, ie this is \code{CNOlist$valueSignals}
}
  \item{times}{
a vector of times, its length should be the same as the number of fields in SimResults and ExpResults
}
  \item{namesCues}{
a vector of names, typically \code{CNOlist$namesCues}
}
  \item{namesSignals}{
a vector of names, typically \code{CNOlist$namesSignals}
}
  \item{valueCues}{
a matrix of dimensions (number of conditions) * (number of cues), typically \code{CNOlist$valueCues}
}
  \item{filename}{
a name for your file, eg. "plot.pdf"
}
}
\details{
The coloring of the background is done as follows: the mean absolute difference between observed and simulated values are computed, and colours are chosen based on this value: red (above 0.9), indianred1 (between O.8 and 0.9), lightpink2 (between 0.7 and 0.8), lightpink (between 0.6 and 0.7), mistyrose (between 0.5 and 0.6), palegoldenrod (between 0.4 and 0.5), palegreen (between 0.3 and 0.4), darkolivegreen3 (between 0.2 and 0.3), chartreuse3 (between 0.1 and 0.2), forestgreen (between 0 and 0.1). This function is used inside \code{cutAndPlotResultsT1}.
}
\value{
This function doesn't return anything, it just produces a plot in a pdf document in your working directory.
}
\author{
C. Terfve
}
\seealso{
plotOptimResults, cutAndPlotResultsT1
}
\examples{
tmpdir<-tempdir()
setwd(tmpdir)

#We will plot the fit of the full initial model compared to the data, without any optimisation
#This is normally not done on a stand alone basis, but if you have a model and would like to visualise 
#its output compared to your data, then this is what you should do

#load and prepare data

data(CNOlistToy,package="CellNOptR")
data(ToyModel,package="CellNOptR")
indicesToy<-indexFinder(CNOlistToy,ToyModel,verbose=TRUE)
ToyFields4Sim<-prep4Sim(ToyModel)

#simulate the model

Sim<-simulatorT1(CNOlist=CNOlistToy,Model=ToyModel,SimList=ToyFields4Sim,indexList=indicesToy)

#format the results and data as expected by plotOptimResults

SimRes<-Sim[,indicesToy$signals]
SimResults<-list(t0=matrix(data=0,nrow=dim(SimRes)[1],ncol=dim(SimRes)[2]),t1=SimRes)
expResults<-list(t0=CNOlistToy$valueSignals[[1]],t1=CNOlistToy$valueSignals[[2]])

#plot

plotOptimResultsPDF(
	SimResults=SimResults,
	expResults=expResults,
	times=CNOlistToy$timeSignals[1:2],
	namesCues=CNOlistToy$namesCues,
	namesSignals=CNOlistToy$namesSignals,
	valueCues=CNOlistToy$valueCues, 
	filename="Toyfull.pdf")
}

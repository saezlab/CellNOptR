\name{plot-method}
\docType{methods}
\alias{plot.CNOlist}
\alias{plot,CNOlist,ANY-method}



\title{plot a \code{"CNOlist"} object - methods}
\description{
 A plot method for \code{CNOlist}. 
}

\section{Usage}{
    \code{signature(x="CNOlist")}: Please see the help page for the \code{plot.CNOlist} method in the \code{CellNOptR} package
}
\section{arguments}{
    \describe{
	\item{x}{The \code{CNOlist} object to plot}
    }
}


\author{
T.Cokelaer
}
\seealso{
\code{\link{readMIDAS}}, \code{\link{makeCNOlist}}
}
\examples{
showClass("CNOlist")

files<-dir(system.file("ToyModel",package="CellNOptR"),full=TRUE)
cnolist  = CNOlist(files[[1]])
# accessors:
cnolist@cues
cnolist@inhibitors
cnolist@signals
cnolist@timepoints
cnolist@stimuli
# In version 1.3.30 and above, use the plot method instead of former plotCNOlist function.
plot(cnolist)
}
\keyword{methods}
\keyword{CNOlist}

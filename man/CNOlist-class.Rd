\name{CNOlist-class}
\docType{class}
\alias{CNOlist-class}
\alias{CNOlist}
\alias{show,CNOlist-method}
\alias{plot,CNOlist-method}

\title{Class \code{"CNOlist"}}
\description{
This function takes as input the filename of a MIDAS file (or the list returned
by \link{makeCNOlist}) and returns an instance of CNOlist class. It provides an 
object oriented approach to manipulate CNOlist. This function calls
\link{readMIDAS} and \link{makeCNOlist}.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("CNOlist", ...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{cues}:}{Object of class \code{"matrix"} }
    \item{\code{inhibitors}:}{Object of class \code{"matrix"}  }
    \item{\code{stimuli}:}{Object of class \code{"matrix"} }
    \item{\code{signals}:}{Object of class \code{"list"}  }
    \item{\code{timepoints}:}{Object of class \code{"vector"} timepoints contained in the signals matrix. }
  }
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object="CNOlist")}: ... }
%    \item{plot}{\code{signature(x="CNOlist", y="missing")}: ... }
    \item{plot}{\code{signature(x="CNOlist")}: Please see the help page for the \code{plot.CNOlist}  method in the \code{CellNOptR} package}

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
\keyword{classes}

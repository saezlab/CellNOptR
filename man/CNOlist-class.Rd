\name{CNOlist-class}
\Rdversion{1.1}
\docType{class}
\alias{CNOlist-class}
\alias{CNOlist}

\title{Class \code{CNOlist}}
\description{
This function takes as input the filename of a MIDAS file and returns an
instance of CNOlist class. It provides an object oriented approach to manipulate
CNOlist. This function calls \link{readMIDAS} and \link{makeCNOlist}.
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
No methods defined with class "CNOlist" in the signature.
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

}
\keyword{classes}

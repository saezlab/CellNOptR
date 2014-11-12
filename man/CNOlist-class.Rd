\name{CNOlist-class}
\docType{class}
\alias{CNOlist-class}
%\alias{CNOlist,ANY}
\alias{CNOlist}


\alias{compatCNOlist,CNOlist-method}
\alias{getCues,CNOlist-method}
\alias{getInhibitors,CNOlist-method}
\alias{getSignals,CNOlist-method}
\alias{getStimuli,CNOlist-method}
\alias{getTimepoints,CNOlist-method}
\alias{getVariances,CNOlist-method}
\alias{plot,CNOlist,CNOlist-method}
%\alias{plot,CNOlist,missing-method}
%\alias{plot,CNOlist,CNOlist,CNOlist-method}
\alias{length,CNOlist-method}
\alias{randomize,CNOlist-method}
\alias{setSignals<-,CNOlist-method}


%\alias{plot,CNOlist-method}
%\alias{plot,CNOlist,ANY,CNOlist-method}
%\alias{plot,CNOlist,CNOlist,CNOlist-method}
%\alias{plot,CNOlist,ANY}
%\alias{plot,CNOlist,CNOlist}







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
    \item{\code{variances}:}{Object of class \code{"list"}  }
    \item{\code{timepoints}:}{Object of class \code{"vector"} timepoints
    contained in the signals matrix. }

  See \link{CNOlist-methods} for details
}
}


\section{Methods}{
  Available methods are plot, compatCNOlist, randomize, length.  See
\link{CNOlist-methods} for details.
}


\author{
T. Cokelaer
}
\seealso{
\code{\link{randomizeCNOlist}}, \code{\link{plotCNOlist}}, \code{\link{plotCNOlist2}}
}
\examples{
showClass("CNOlist")

files<-dir(system.file("ToyModel",package="CellNOptR"),full=TRUE)
cnolist  = CNOlist(files[[1]])
# getters:
getCues(cnolist)
getInhibitors(cnolist)
getSignals(cnolist)
getVariances(cnolist)
getTimepoints(cnolist)
getStimuli(cnolist)
# In version 1.3.30 and above, use the plot method instead of former plotCNOlist function.
plot(cnolist)
new_cnolist = randomize(cnolist)
length(cnolist)
}
\keyword{classes}

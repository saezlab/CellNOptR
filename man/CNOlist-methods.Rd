\name{CNOlist-methods}
\docType{methods}
\alias{CNOlist-methods}

\alias{getCues}
\alias{getCues,CNOlist}

\alias{getSignals}
\alias{getSignals,CNOlist}
%\alias{getSignals,CNOlist,CNOlist-method}

\alias{compatCNOlist}
\alias{compatCNOlist,CNOlist}
%\alias{compatCNOlist,CNOlist,CNOlist-method}

\alias{getStimuli}
\alias{getStimuli,CNOlist}
%\alias{getStimuli,CNOlist,CNOlist-method}

\alias{getTimepoints}
\alias{getTimepoints,CNOlist}
%\alias{getTimepoints,CNOlist,CNOlist-method}

\alias{getVariances}
\alias{getVariances,CNOlist}
\alias{getVariances,CNOlist,CNOlist-method}

\alias{setSignals<-}
\alias{setSignals<-,CNOlist}
\alias{setSignals<-,CNOlist,CNOlist-method}

\alias{getInhibitors}
\alias{getInhibitors,CNOlist}
%\alias{getInhibitors,CNOlist-method}

\alias{randomize}
\alias{randomize,CNOlist}
\alias{randomize,CNOlist,CNOlist-method}

\alias{length}
\alias{length,CNOlist}
\alias{length,CNOlist,ANY-method}

%\alias{plot}
%\alias{plot,CNOlist}
\alias{plot,CNOlist,CNOlist,CNOlist-method}
%\alias{plot,CNOlist,CNOlist,CNOlist-methods}
%\alias{plot,CNOlist,missing,CNOlist-method}



\title{List of \link{CNOlist-class} methods}
\description{
 CNOlist is a class with a set of methods described here below.

}

\section{Usage}{
    \code{signature(x="CNOlist")}
}
\section{Getters}{
    \describe{
	\item{getCues}{Returns the cues (matrix) found in the \code{CNOlist}}
	\item{getSignals}{Returns the signals (list of matrices) found in the \code{CNOlist}}
	\item{getStimuli}{Returns the cues found in the \code{CNOlist}}
	\item{getInhibitors}{Returns the inhibitors found in the \code{CNOlist}}
	\item{getTimepoints}{Returns the timepoints found in the \code{CNOlist}}
	\item{getVariances}{Returns the Variances (list of matrices) found in the
        \code{CNOlist}. Will be different from zero only if replicates were found in the
        MIDAS data. See \code{makeCNOlist}}
    }
}
\section{Setters}{
    \describe{
    \item{setSignals}{Set the signals. No sanity check!}
    }
}

\section{Other methods}{
    \describe{
      \item{compatCNOlist}{convert the instance CNOlist into the old-style
            returned by \code{\link{makeCNOlist}} that is a list. Used in the ODE package. }

      \item{length}{returns length of CNOlist (number of time points)}

      \item{randomize}{randomizes the signals matrice in a CNOlist. See
            \code{\link{randomizeCNOlist}} for details}

      \item{show}{prints summary information }

      \item{plot}{plot the CNOlist instance using the \link{plotCNOlist} function.}

      \item{plot}{\code{signature(x="CNOlist", y="CNOlist")}: Please see the
          page of \code{\link{plotCNOlist2}} for more details.}
      }
}



\author{
T.Cokelaer
}
\seealso{
\code{\link{CNOlist-class}}, \code{\link{randomizeCNOlist}} \code{\link{makeCNOlist}}
}
\examples{
showClass("CNOlist")

data(CNOlistToyMMB, package="CellNOptR")
cnolist = CNOlistToyMMB

# In version 1.3.30 and above, use the plot method instead of former plotCNOlist function.
plot(cnolist)
# In version 1.5.14 and above, use getters instead of the @ operator
getCues(cnolist)

# others:
new_cnolist = randomize(cnolist)

}
\keyword{methods}
\keyword{CNOlist}

\name{checkSignals}
\alias{checkSignals}
\title{
Check the CNOlist and model matching
}
\description{
This function checks that all the signals in a CNOlist match to species in the 
model. It also checks that the CNOlist and Model lists have the right format and 
contain the right fields. It is called by the \link{preprocessing} function so there is 
no need to call it directly anymore if you use the \link{preprocessing} function.

In version 1.3.20, check of inhibitors and stimuli is also performed. 
}
\usage{
checkSignals(CNOlist, model)
}
\arguments{
    \item{CNOlist}{
        A CNOlist structure, as created by \link{makeCNOlist}.
    }
    \item{model}{
        A model structure, as created by \link{readSIF}.
    }
}
\details{
If the formats are wrong, this function produces an error. If the
signals/inhibitors/stimuli do not match the species, this function produces a
warning that explains which signals does no match any species.
}

\author{
C. Terfve, T. Cokelaer
}
\seealso{
\code{\link{makeCNOlist}}, \link{readSIF}, \link{preprocessing}
}
\examples{
data(CNOlistToy, package="CellNOptR")
data(ToyModel, package="CellNOptR")
checkSignals(CNOlistToy, ToyModel)
}



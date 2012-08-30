\name{preprocessing}
\alias{preprocessing}
\title{
Performs the pre-processing steps
}
\description{
This function performs any of the following preprocessing steps:
    \enumerate{
        \item{removes Non-Controllable and Non-Observables nodes}
        \item{compress the model}
        \item{and-gates expansion}
    }
}
\usage{
preprocessing(data, model, cutNONC=TRUE, compression=TRUE,
    expansion=TRUE, ignoreList=NA, maxInputsPerGate=2, verbose=TRUE)
}

\arguments{
  \item{data}{
the CNOlist that contains the data that you will use
}
  \item{model}{
the model object as returned by \link{readSIF}
}
  \item{cutNONC}{
Removes the NONC nodes using \link{cutNONC} and \link{findNONC} (Default is TRUE).
}
  \item{compression}{
Compress the model using \link{compressModel} (Default is TRUE).
}
  \item{expansion}{
Add and gates using \link{expandGates} (Default is TRUE).
}
  \item{ignoreList}{
list of reactions indices to ignore in the expansion step.
}
  \item{maxInputsPerGate}{
used by the \link{expandGates} function to set maximum inputs per and gates.
}
  \item{verbose}{
verbose option (Default is TRUE).
}
}


\details{
  The function can apply any or none of the pre-processing steps. It returns the
new model and the indices returned by indexFinder.

}

\value{
  the new model
}

\author{
 T. Cokelaer
}


\seealso{
\link{readSIF}, \link{readMIDAS}, \link{cutNONC}, \link{findNONC},
\link{compressModel}, \link{expandGates}.

}
\examples{
data(ToyModel,package="CellNOptR")
data(CNOlistToy,package="CellNOptR")
model = preprocessing(CNOlistToy, ToyModel, cutNONC=FALSE)
}

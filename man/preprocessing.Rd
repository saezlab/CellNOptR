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
preprocessing(Data, Model, cutnonc=TRUE, compression=TRUE,
expansion=TRUE,ignoreList=NA, maxInputsPerGate=2, verbose=TRUE)
}

\arguments{
  \item{Data}{
   the CNOlist that contains the data that you will use
  }
  \item{Model}{
    the model object as returned by \link{readSif}
  }
  \item{cutnonc}{
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
  \item{model}{the new model}
  \item{indices}{the indices of the new model found in the CNOlist object using
indexFinder.}
}

\author{
 T. Cokelaer
}


\seealso{
\link{readSif}, \link{readMIDAS}, \link{cutNONC}, \link{findNONC},
\link{compressModel} and \link{expandGates}.

}
\examples{
data(ToyModel,package="CellNOptR")
data(CNOlistToy,package="CellNOptR")
res = preprocessing(CNOlistToy, ToyModel, cutnonc=FALSE)
# res is a list containing the new model (res$model) and its indices
# (res$indices).
}

\name{cutCNOlist}
\alias{cutCNOlist}
\title{
Cut a CNOlist structure according to a model
}
\description{
The MIDAS file may contain species that are not contained in the model. If you
want to remove cues and signals from your CNOlist that are not contained in the
model, you can use this function by using the parameter model. 
It can also be used to remove some time points by using the parameter
cutTimeIndices. Both parameters can be used at the same time and at least one of
them must be provided.
}
\usage{
    cutCNOlist(cnolist, model, cutTimeIndices, verbose=FALSE)
}

\arguments{
  \item{cnolist}{the CNOlist structure}
  \item{model}{the model}
  \item{cutTimeIndices}{the time indices to remove}
  \item{verbose}{Set it to True to get the signals and cues not found in the model}
}


\value{
  \item{cutCNOlist}{the new CNOlist object}
}

\seealso{
    \link{CNOlist-class}
}

\note{added in version 1.5.10}

\author{
 T. Cokelaer
}



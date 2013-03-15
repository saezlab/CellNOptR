\name{cutCNOlist}
\alias{cutCNOlist}
\title{
Cut a CNOlist structure according to a model
}
\description{
The MIDAS file may contain species that are not contained in the model. If you
want to remove cues and signals from your CNOlist that are not contained in the
model, you can use this function. 
}
\usage{
    cutCNOlist(cnolist, model, verbose=FALSE)
}

\arguments{
  \item{cnolist}{the CNOlist structure}
  \item{model}{the model}
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



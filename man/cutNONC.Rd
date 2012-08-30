\name{cutNONC}
\alias{cutNONC}
\title{
Cuts the non-observable/non-controllable species from the model
}
\description{
This function cuts the non-observable and/or non-controllable species from the model, and returns a cut model.
}
\usage{
cutNONC(model, NONCindexes)
}
\arguments{
  \item{model}{
a model structure, as produced by readSIF
}
  \item{NONCindexes}{
a vector of indices of species to remove in that model, as produced for example by findNONC
}
}
\details{
This function takes in a model and a vector of indices of species to remove in that model and it removes those species and any reaction involving them (be aware, if you have x&y=z and x is to be removed, then the function produces y=z, because it works by removing entire rows of the model matrices and then removes the columns that do not have either an input or an output). This function could actually be used to cut any species, not only NONC species.
}
\value{
a model
}
\author{
C.Terfve
}
\seealso{
\link{findNONC}, \link{readSIF}
}
\note{No need to call this function directly since version 0.99.24. Use
\link{preprocessing} instead.}

\examples{
data(CNOlistToy,package="CellNOptR")
data(ToyModel,package="CellNOptR")
indicesToy<-indexFinder(CNOlistToy,ToyModel,verbose=FALSE)
ToyNCNOindices<-findNONC(ToyModel,indicesToy,verbose=FALSE)
ToyNCNOcut<-cutNONC(ToyModel,ToyNCNOindices)
}

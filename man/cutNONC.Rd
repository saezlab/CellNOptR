\name{cutNONC}
\alias{cutNONC}
\title{
Cuts the non-observable/non-controllable species off the model
}
\description{
This function cuts the non observable and/or non controllable species off the model, and returns a cut model
}
\usage{
cutNONC(Model, NONCindexes)
}
\arguments{
  \item{Model}{
a model structure, as produced by readSif
}
  \item{NONCindexes}{
a vector of indices of species to remove in that model, as produced for example by findNONC
}
}
\details{
This function takes in a model and a vector of indices of species to remove in that model and it removes those species and any reaction involving them (be aware, if you have x&y=z and x is to be removed, then the function produces y=z, because it works by removing entire rows of the Model matrices and then removes the columns that do not have either an input or an output). This function could actually be used to cut any species, not only NONC species.
}
\value{
a model
}
\author{
C.Terfve
}
\seealso{
findNONC, readSif
}
\examples{
data(CNOlistToy,package="CellNOptR")
data(ToyModel,package="CellNOptR")
indicesToy<-indexFinder(CNOlistToy,ToyModel,verbose=FALSE)
ToyNCNOindices<-findNONC(ToyModel,indicesToy,verbose=FALSE)
ToyNCNOcut<-cutNONC(ToyModel,ToyNCNOindices)
}

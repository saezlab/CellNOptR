\name{findNONC}
\alias{findNONC}
\title{
Find the indexes of the non-observable and non controllable species
}
\description{
This function finds the indexes of the non-observable and non controllable species and returns the indices, in the model, of the species to remove
}
\usage{
findNONC(Model, indexes, verbose)
}
\arguments{
  \item{Model}{
a model, as created by readSif
}
  \item{indexes}{
 a list of indices of species measured, stimulated or inhibited, as created by indexFinder
}
  \item{verbose}{
 do you want information about the ncno species printed on the screen? Default to FALSE but we would advise to put it to true the first time that the function is called
}
}
\details{
This function uses the function floyd.warshall.all.pairs.sp from the package RBGL.  Non observable nodes are those that do not have a path to any measured species in the model, whereas non controllable nodes are those that do not receive any information from a species that is perturbed in the data.
}
\value{
a vector of indices of species to remove
}
\author{
C. Terfve
}
\seealso{
cutNONC, indexFinder, readSif
}
\examples{
data(CNOlistToy,package="CellNOptR")
data(ToyModel,package="CellNOptR")
checkSignals(CNOlistToy,ToyModel)
indicesToy<-indexFinder(CNOlistToy,ToyModel,verbose=FALSE)
ToyNCNOindices<-findNONC(ToyModel,indicesToy,verbose=FALSE)
}

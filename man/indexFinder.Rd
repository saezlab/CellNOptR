\name{indexFinder}
\alias{indexFinder}
\title{
Finds the indices, in the Model fields, of the species that are measured/inhibited/stimulated
}
\description{
This function finds the indices, in the Model fields, of the species that are measured/inhibited/stimulated. It looks for their position in Model$namesSpecies which has the same order as the rows of interMat and notMat, and therefore these indexes can be used there as well.
}
\usage{
indexFinder(CNOlist, Model,verbose=FALSE)
}
\arguments{
  \item{CNOlist}{
a CNOlist structure, as produced by makeCNOlist
}
  \item{Model}{
a model structure, as produced by readSif
}
  \item{verbose}{
do you want information about the cues and signals identities printed on the screen? Default if false but we would advise to set it to true when the function is called for the first time.
}
}
\value{
a list with fields:
\item{signals}{vector of indices of the measured species}
\item{stimulated}{vector of indices of the stimulated species}
\item{inhibited}{vector of indices of the inhibited species}
}
\author{
C. Terfve
}
\seealso{
makeCNOlist, ReadSif
}
\examples{
data(CNOlistToy,package="CellNOptR")
data(ToyModel,package="CellNOptR")
indicesToy<-indexFinder(CNOlistToy,ToyModel,verbose=TRUE)
}

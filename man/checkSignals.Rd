\name{checkSignals}
\alias{checkSignals}
\title{
Check the CNOlist and model matching
}
\description{
This function checks that all the signals in the CNOlist match to species in the model. It also checks that the CNOlist and Model lists have the right format and contain the right fields.
}
\usage{
checkSignals(CNOlist, Model)
}
\arguments{
  \item{CNOlist}{
A CNOlist structure, as created by makeCNOlist
}
  \item{Model}{
A model structure, as created by readSif
}
}
\details{
If the formats are wrong, this function produces an error with an explanation message. If the signals don't match the species, this function produces a warning that explains which signals don't match any species, and advises to remove them (THIS SHOULD BE DONE IMPERATIVELY) which is not done at this point but could be done in the form of this function returning a new CNOlist which would be the original if all ok or a CNOlist with non matching signals removed when necessary. I don't see this as a priority at this stage so if it happens, the user should remove those signals manually from the CNOlist (in both fields valueSignals and namesSignals).
}
\value{
If all ok this function doesn't return anything.
}
\author{
C. Terfve
}
\seealso{
makeCNOlist, readSif
}
\examples{
data(CNOlistToy,package="CellNOptR")
data(ToyModel,package="CellNOptR")
checkSignals(CNOlistToy,ToyModel)
}



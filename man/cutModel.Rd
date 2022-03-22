\name{cutModel}
\alias{cutModel}
\title{
Cut a model structure according to a bitstring
}
\description{
This function is for developers only. 
}
\usage{
cutModel(model, bString)
}

\arguments{
  \item{model}{the model object to cut.}
  \item{bString}{the bitString to be used to cut the model object.}
}


\value{
  \item{cutModel}{the new model object}
}

\note{added in version 1.3.16}

\author{
 T. Cokelaer
}

\examples{
   
    data(ToyModel,package="CellNOptR")
    
    bString = rep(1,length(ToyModel$reacID))
    # remove some reactions by setting them to 0. 
    bString[c(2,5,6)] <- 0
    
    prepModel = cutModel(model = ToyModel, bString = bString )
 
    
}
\name{defaultParameters}
\alias{defaultParameters}
\title{
Create a list of default parameters
}
\description{
This function provides a list of default parameters including the Genetic
Algorithm parameters.}
\usage{
defaultParameters(data=NA, model=NA)
}

\arguments{
  \item{data}{
   a CNOlist structure, as created by \code{makeCNOlist}
  }
  \item{model}{
    a model structure, as created by \code{readSIF}, normally pre-processed but 
    that is not a requirement of this function
  }
}


\details{
  The list contains the Genetic Algorithm parameter, a verbose option and can be
  used to store the Data and Model.
}

\value{
  \item{params}{a list with the fields: \code{data}, \code{model},
\code{verbose} and all default parameters of \link{gaBinaryT1} }
}
\author{
 T. Cokelaer
}
\examples{
data(ToyModel, package="CellNOptR")
data(CNOlistToy, package="CellNOptR")
params = defaultParameters(CNOlistToy, ToyModel)
}

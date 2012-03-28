\name{defaultParameters}
\alias{defaultParameters}
\title{
Create a list of default parameters
}
\description{
This function provides a list of default parameters including the Genetic
Algorithm parameters.}
\usage{
defaultParameters(Data=NA, Model=NA)
}

\arguments{
  \item{Data}{
   the CNOlist that contains the data that you will use
  }
  \item{Model}{
    the model that you want to optimise
  }
}


\details{
  The list contains the Genetic Algorithm parameter, a verbose option and can be
  used to store the Data and Model.
}

\value{
  \item{params}{a list of default parameters. }
}
\author{
 T. Cokelaer
}
\examples{
data(ToyModel,package="CellNOptR")
data(CNOlistToy,package="CellNOptR")
params = defaultParameters(CNOlistToy, ToyModel)
}

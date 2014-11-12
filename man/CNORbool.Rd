\name{CNORbool}
\alias{CNORbool}
\title{
Simple Boolean analysis standalone
}
\description{
This function performs the optimisation of a PKN model to a CNOlist data set. 
It optimises each time point found in the data and returns the processed model
as well as a list of optimised bitstring corresponding to each time points that
has been optimised.

This function does not create any plots or reports unlike \code{\link{CNORwrap}}.
}
\usage{
CNORbool(CNOlist, model, paramsList=defaultParameters(),
    compression=TRUE, expansion=TRUE, cutNONC=TRUE, verbose=FALSE,
    timeIndices = NULL)
}

\arguments{
\item{CNOlist}{
a CNOlist structure, as created by \link{makeCNOlist} or a MIDAS filename
}
  \item{model}{
a model structure, as created by \link{readSIF} or a SIF filename.
}

  \item{paramsList}{
   Parameters of the genetic algorithm. If not provided, it is populated with the \link{defaultParameters} function.
   }
\item{compression}{
compress the model (default TRUE)
}
\item{expansion}{
expand the gates (default TRUE)

}
\item{cutNONC}{
cut the NONC nodes off the model by (default TRUE)
}
\item{verbose}{
FALSE}
\item{timeIndices}{by default, optimise T1 and T2 assuming there are the 2 first
time points. However, with this argument you can change that behaviour to
arbitrary time points.}

}
\value{
This function returns 2 components. The first one is the processed model used in
the optimisation. The second is a list of optimised bitstrings found for each
time points available in the MIDAS data set.
}
\author{
T.Cokelaer, S.Schrier
}
\examples{

data(CNOlistToy,package="CellNOptR")
data(ToyModel,package="CellNOptR")
res = CNORbool(CNOlist=CNOlistToy, model=ToyModel)

}

\name{readSBMLQual}
\alias{readSBMLQual}
\alias{readSBMLQual}
\title{
Read a SBMLQual document and returns a SIF object (as returned by readSIG
}
\description{
This function reads a SBMLQual XML file where a model is stored. The XML is
scanned and saved in SIF format in a temporary file. This file is read by
readSIF. The returned object is therefore the output of readSIF.
}
\usage{
readSBMLQual(filename)
}
\arguments{
  \item{filename}{
The name of a SBMLQual file.
}
}
\value{
a model list with fields:
    \item{interMat}{contains a matrix with column for each reaction and a row for each species, with a -1 where the species is the source node and a +1 where the species is a target node, and 0 otherwise}
    \item{notMat}{has the same format as interMat but just contains a 1 if the source node enters the reac with a negative effect, and 0 otherwise}
    \item{namesSpecies}{vector that contains the names of the species in the same order as the rows of the interMat and notMat matrices}
    \item{reacID}{vector that holds character strings specifying the reaction in full letters, in the same order as the columns of interMat and notMat}
}

\references{
SBMLQual: qualitative models. See SBML.org for details.
}
\author{
T. Cokelaer
}
\examples{

\dontrun{
sif = readSBMLQual("test.xml")
}
}

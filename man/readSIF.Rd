\name{readSIF}
\alias{readSif}
\alias{readSIF}
\title{
Read a SIF file and create a model object
}
\description{
This function reads in a cytoscape SIF file and creates a model object that can be used in the CellNOptR procedure.
}
\usage{
readSIF(sifFile)
}
\arguments{
  \item{sifFile}{
The name of a SIF file (Cytoscape format). See details for the accepted format.
}
}
\details{
The input argument is the name of a SIF file that contains a prior knowledge
network. The SIF format (See http://wiki.cytoscape.org/Cytoscape_User_Manual/Network_Formats/) 
can be of the form:

\code{ nodeA typeA node2}

or

\code{ node2 typeB node3 node4}

with space or tabulations in between types and nodes. Spaces and tabulations
have different meaning in the original cytoscape format. However, we do not
differentiate them. Therefore names of the nodes cannot have white space inside
them.

The accepted values for the types are -1 (inhibitor) or 1 (normal relation).

If there are ANDs they should be introduced as dummy nodes called "and#" (don't
forget the number after "and" otherwise this won't be recognised).  Please be
aware that "and" nodes are not expected to be negated, i.e. there are not
supposed to be !and1=xyz because that amounts to inverting the sign of all
inputs of and1, which is more simply done at the inputs level.
}
\value{
a model list with fields:
    \item{interMat}{contains a matrix with column for each reaction and a row for each species, with a -1 where the species is the source node and a +1 where the species is a target node, and 0 otherwise}
    \item{notMat}{has the same format as interMat but just contains a 1 if the source node enters the reac with a negative effect, and 0 otherwise}
    \item{namesSpecies}{vector that contains the names of the species in the same order as the rows of the interMat and notMat matrices}
    \item{reacID}{vector that holds character strings specifying the reaction in full letters, in the same order as the columns of interMat and notMat}
}

\references{
Shannon P, Markiel A, Ozier O, Baliga NS, Wang JT, Ramage D, Amin N, Schwikowski
B, Ideker T. Cytoscape: a software environment for integrated models of
biomolecular interaction networks. Genome Research 2003 Nov; 13(11):2498-504.
}
\author{
C. Terfve
}
\examples{
tmpdir<-tempdir()
setwd(tmpdir)
cpfile<-dir(system.file("ToyModel",package="CellNOptR"),full=TRUE)
file.copy(from=cpfile,to=getwd(),overwrite=TRUE)
ToyModel<-readSIF(sifFile="ToyPKNMMB.sif")
}

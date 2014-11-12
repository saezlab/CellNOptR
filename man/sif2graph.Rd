\name{sif2graph}
\alias{sif2graph}

\title{
  Convert sif to graph
}

\description{
  This function receives as input a network in form sif format and converts it to graph format.
}

\usage{
sif2graph(sif)
}


\arguments{
  \item{sif}{the name of a sif file or the equivalent table}
}

\details{
This function takes a network in sif format (tabel or file), i.e. sourceNode-tab-sign-tab-targetNode. If there are ANDs they should be introduced as dummy nodes called and# (don't forget the number after "and" otherwise this won't be recognised).  Please be aware that "and" nodes are not expected to be negated, i.e. there are not supposed to be !and1=xyz because that amounts to inverting the sign of all inputs of and1, which is more simply done at the inputs level.

In the resulting graph, the sign of each link is encodes as the weigth of the edge (-1 negative regulation, +1 positive regulation).
}

\value{
\item{g}{the corresponding graph}
}

\author{
    F.Eduati
}

\seealso{
\link{graph2sif}, \link{model2sif}, \link{readSIF}, 
}



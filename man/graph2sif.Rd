\name{graph2sif}
\alias{graph2sif}

\title{
Convert graph to SIF
}

\description{
  This function converts a network form graph format to SIF format. The resulting table can also be saved in a SIF file.
}

\usage{
graph2sif(graph, writeSif=FALSE, filename="Graph")

}


\arguments{
  \item{graph}{
a graph, as generated using \link{sif2graph}
}
  \item{writeSif}{
if writeSif=FALSE (default) the SIF file is not saved. If writeSif=TRUE it is saved.
}
  \item{filename}{
the name of the SIF file saved if writeSif=TRUE. Default is Model.sif.
}

}


\details{
The sign of link is supposed to be encoded in the graph as the weigth of the edge (-1 negative regulation, +1 positive regulation).
}

\value{
\item{sifFile}{a table with all the links in the model in the format sourceNode-tab-sign-tab-targetNode}
}

\author{
    F. Eduati
}

\seealso{
\link{model2sif}, \link{sif2graph}, \link{readSIF}, 
}



\name{model2igraph}
\alias{model2igraph}

\title{
  Convert a model object to a igraph object
}

\description{
  This function receives as input a model object and converts it to a graph
object made by igraph. igraph provides lots of utilities especially to write the
file in different format such as GML.
}

\usage{
model2igraph(model)

}


\arguments{
  \item{model}{
the model as generated using \link{readSIF}
}
}



\value{
\item{g}{a igraph object}
}

\author{
    T. Cokelaer
}

\seealso{
\link{graph2sif}, \link{sif2graph} 
}




\name{mapBack}
\alias{mapBack}

\title{
  Map an optimised model back onto the PKN model.
}

\description{
  Map an optimised model back onto the PKN model.
}

\usage{
mapBack(model, PKN, bString)

}


\arguments{
  \item{model}{
the optimised model
}
  \item{PKN}{
the Prior Knowledge network
}
  \item{bString}{
the optimised bitString
}
}



\value{
\item{bStringPKN}{the corresponding bitstring corresponding to the original PKN.
}
}

\author{
    F.Eduati
}

\seealso{
\link{graph2sif}, \link{sif2graph}, \link{readSIF}, 
}



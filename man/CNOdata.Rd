\name{CNOdata}
\alias{CNOdata}
\title{
Get data from a CellNOpt data repository
}
\description{
This function fetch a file from the URL provided (default is
http://www.ebi.ac.uk~/cokelaer/cellnopt/data/_downloads)
and save it into a temporary file.

You will need Rcurl package to use this function.
}
\usage{
CNOdata(filename, verbose=FALSE, url=NULL)
}
\arguments{
  \item{filename}{a valid filename that can be found in the url}
  \item{verbose}{ FALSE by default, it prints the path of the temporary file
where data has been copied.}
   \item{url}{ You can overwrite the default URL
(http://www.ebi.ac.uk~/cokelaer/cellnopt/data/_downloads) with this argument.
}
}
\value{
    the path of the temporary file where data has been copied.
}

\author{
T. Cokelaer
}
\examples{
readSIF(CNOdata("PKN-ToyMMB.sif"))

}



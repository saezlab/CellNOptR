\name{writeMIDAS}
\alias{writeMIDAS}
\title{Write a CNOlist structure into a MIDAS file}
\description{
This function takes a CNOlist structure (output of makeCNOlist and readMIDAS) and save it a MIDAS format.
}
\usage{
writeMIDAS(CNOlist, filename)
}

\arguments{
    \item{CNOlist}{
        a CNOlist structure
    }
    \item{filename}{
        a filename. Not overwritten if it exists already    
    }
}

\author{
T. Cokelaer
}
\seealso{
\link{makeCNOlist}, \link{readMIDAS}, \link{CNOlist}
}
\examples{
   tmpdir<-tempdir()
   setwd(tmpdir)
   data(CNOlistToy)
   writeMIDAS(CNOlistToy, 'test.csv')
   readMIDAS('test.csv')
}

\name{writeMIDAS}
\alias{writeMIDAS}
\title{Write a CNOlist structure into a MIDAS file}
\description{
This function takes a CNOlist structure (output of makeCNOlist and readMIDAS) and save it a MIDAS format.
}
\usage{
writeMIDAS(CNOlist, filename,  timeIndices=NULL, overwrite=FALSE)
}

\arguments{
    \item{CNOlist}{
        a CNOlist structure
    }
    \item{filename}{
        a filename. Not overwritten if it exists already    
    }
    \item{timeIndices}{
        select subset of the times to be saved. Works with indices (not time values)
    }
    \item{overwrite}{
        overwrite the file if it exists already (Default is FALSE)
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
   writeMIDAS(CNOlistToy, 'test.csv', timeIndices=c(1,2), overwrite=TRUE)
}

\name{readMIDAS}
\alias{readMIDAS}
\alias{readMidas}
\title{Reads in a CSV MIDAS file}
\description{
This function takes in a single argument, the name of a CSV MIDAS file
containing the data, and returns a list that contains all the elements to build
a CNOlist. The output of this function should be used as input for
\code{\link{makeCNOlist}}.
}
\usage{
readMIDAS(MIDASfile, verbose=TRUE)
}

\arguments{
    \item{MIDASfile}{
        a CSV MIDAS file (see details)
    }
    \item{verbose}{
        logical (default to TRUE).
    }
}

\details{
This function does not return a CNOlist, but the output of this function can be 
used directly into makeCNOlist to create one. The MIDAS file format is described
in Saez-Rodriguez et al. (2008).

If you have all of the readouts measured at the same series of time points, you
can specify a unique DA: column which must have the format "DA:ALL". 

}
\value{
this function returns a list with fields:
\item{dataMatrix}{matrix containing the data in the MIDAS file}
\item{TRcol}{indexes of the columns that contain the treatments (excluding cell line)}
\item{DAcol}{indexes of the columns that contain the data time points}
\item{DVcol}{indexes of the columns that contain the actual values (measurements)}
}
\references{
J. Saez-Rodriguez, A. Goldsipe, J. Muhlich, L. Alexopoulos, B. Millard, D. A. Lauffenburger, P. K. Sorger 
\emph{Flexible Informatics for Linking Experimental Data to Mathematical Models via DataRail.} Bioinformatics, 24:6, 840-847 (2008). 
}
\author{
C.Terfve
}
\seealso{
\link{makeCNOlist}
}
\examples{
   
   cpfile<-dir(system.file("ToyModel",package="CellNOptR"),full=TRUE)
   file.copy(from=cpfile,to=getwd(),overwrite=TRUE)
   dataToy<-readMIDAS(MIDASfile='ToyDataMMB.csv')
   CNOlistToy<-makeCNOlist(dataset=dataToy,subfield=FALSE)
}


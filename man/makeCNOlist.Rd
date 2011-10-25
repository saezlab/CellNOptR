\name{makeCNOlist}
\alias{makeCNOlist}
\title{
Make a CNOlist structure
}
\description{
This function takes as input the output of readMIDAS and extracts the elements that are needed in a CNO project
}
\usage{
makeCNOlist(dataset, subfield)
}
\arguments{
  \item{dataset}{
output of readMIDAS
}
  \item{subfield}{
TRUE or FALSE, specifies if the column headers contain subfields or not i.e. if I should look for TR:sthg:sthg or just TR:sthg
}
}
\details{
Be aware that most of the functions in this package, including this one, expect the data to contain measurements at time 0, but these should all be equal to zero according to the normalisation procedure that should be used. Therefore, if you have on time point, the files valueSignals contains two matrices, one for t0 and one for t1.

I there are replicate rows in the MIDAS file (i.e. identical cues and identical time), this function averages the values of the measurements for these replicates.

}
\value{
a CNOlist with fields
\item{namesCues}{a vector of names of cues}
\item{namesStimuli}{a vector of names of stimuli}
\item{namesInhibitors}{a vector of names of inhibitors}
\item{namesSignals}{a vector of names of signals}
\item{timeSignals}{a vector of times}
\item{valueCues}{a matrix of dimensions nConditions x nCues, with 0 or 1 if the cue is present or absent in the particular condition}
\item{valueInhibitors}{a matrix of dimensions nConditions x nInhibitors, with 0 or 1 if the inhibitor is present or absent in the particular condition}
\item{valueStimuli}{of dimensions nConditions x nStimuli, with 0 or 1 if the stimuli is present or absent in the particular condition}
\item{valueSignals}{a list of the same length as timeSignals, each element containing a matrix of dimensions nConditions x nsignals, with the measurements.}
}
\references{
J. Saez-Rodriguez, L. G. Alexopoulos, J. Epperlein, R. Samaga, D. A. Lauffenburger, S. Klamt and P. K. Sorger. Discrete logic modeling as a means to link protein signaling networks with functional analysis of
mammalian signal transduction, Molecular Systems Biology, 5:331, 2009.
}
\author{
C. Terfve
}
\seealso{
readMIDAS
}
\examples{
tmpdir<-tempdir()
setwd(tmpdir)
cpfile<-dir(system.file("ToyModel",package="CellNOptR"),full=TRUE)
file.copy(from=cpfile,to=getwd(),overwrite=TRUE)
dataToy<-readMIDAS(MIDASfile='ToyDataMMB.csv')
CNOlistToy<-makeCNOlist(dataset=dataToy,subfield=FALSE)
}


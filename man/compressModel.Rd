\name{compressModel}
\alias{compressModel}
\title{
Compress a model
}
\description{
This function compresses a model by compressing species that are not 
signals/inhibited/stimulated and that are not dead ends/in complex logic 
(i.e. only species with either one input or one output are compressed)/in self loops.

You can also use \link{preprocessing} function instead that calls \link{compressModel} and 
other preprocessing functions.  
}
\usage{
compressModel(model, indexes)
}
\arguments{
    \item{model}{
        a model structure as produced by \link{readSIF}.
    }
    \item{indexes}{
        list of indexes of the species stimulated/inhibited/measured in the model, as created by \link{indexFinder}.
    }
}
\details{
Be aware that in the multiple inputs/one output case, if one of the outputs is 
an '&' gate this function handles it fine as long as it is an '&' with 2 inputs 
and no more.
}
\value{
a compressed model list, with an additional field called 'speciesCompressed' 
that contains the names of the species that have been compressed
}
\author{
C. Terfve
}
\seealso{
\link{indexFinder}, \link{readSIF}, \link{preprocessing}
}
\note{No need to call this function directly since version 0.99.24. Use 
 \link{preprocessing} instead.}

\examples{
#load data

data(CNOlistToy,package="CellNOptR")
data(ToyModel,package="CellNOptR")

indicesToy<-indexFinder(CNOlistToy,ToyModel,verbose=FALSE)
toyComp<-compressModel(ToyModel,indicesToy)
}

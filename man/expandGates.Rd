\name{expandGates}
\alias{expandGates}
\title{
Expand the gates of a model
}
\description{
This function takes in a model and splits all AND gates into ORs.  In addition,
wherever there are more than one, it creates all possible ANDs combinations of
them, but considering only ANDs with 2, 3 or 4 inputs according to the user
argument (default is 2)
}
\usage{
expandGates(model, ignoreList=NA, maxInputsPerGate=2)
}
\arguments{
  \item{model}{a model structure}
  \item{ignoreList}{a list of reactions indices to ignore}
  \item{maxInputsPerGate}{maximum number of input per gates (Default is 2; up
to 4)}
}
\details{
This function returns a model with additional fields that help keep track of the processing done on the network.  I would advice not to overwrite on the initial model but rather to assign the result of this function to a variable with a different name.
}
\value{
returns a model, with additional fields:
\item{SplitANDs}{list that contains a named element for each AND reac that has been split, and each element contains a vector with the names of the of the reactions that result from the split if nothing was split, this element has the default value $initialReac  [1] "split1" "split2"}
\item{newANDs}{list that contains an element for each new '&' gate, named by the name of this new and reac, and containing a vector of the names of the reactions from which it was created (contains all the reacs in that pool, not the particular ones, this could be improved)}
}
\author{
C.Terfve. T. Cokelaer, A.MacNamara, Martin-Franz-Xaver Pirkl
}
\note{No need to call this function directly since version 0.99.24. Use 
 \link{preprocessing} instead.}


\examples{
#load data

data(CNOlistToy,package="CellNOptR")
data(ToyModel,package="CellNOptR")

#pre-process the model

indicesToy<-indexFinder(CNOlistToy,ToyModel,verbose=TRUE)
ToyNCNOindices<-findNONC(ToyModel,indicesToy,verbose=TRUE)
ToyNCNOcut<-cutNONC(ToyModel,ToyNCNOindices)
indicesToyNCNOcut<-indexFinder(CNOlistToy,ToyNCNOcut)
ToyNCNOcutComp<-compressModel(ToyNCNOcut,indicesToyNCNOcut)
indicesToyNCNOcutComp<-indexFinder(CNOlistToy,ToyNCNOcutComp)
ToyNCNOcutCompExp<-expandGates(ToyNCNOcutComp, maxInputsPerGate=4)
}

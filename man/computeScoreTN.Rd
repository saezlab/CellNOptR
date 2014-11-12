\name{computeScoreTN}
\alias{computeScoreTN}

\title{
  Compute the score at TN of a model/data set using a bitString to cut the model.
}

\description{
The bitString made of 0 and 1 allows to select a submodel from the model provided. 
Then, the simulator function are called to compute the objective function. The 
sizeFac and NAFac are penalties added to the final score as described in 
\code{gaBinaryTN}.
}

\usage{
computeScoreTN(CNOlist, model, simList=NULL, indexList=NULL, simResPrev=NULL,
    bStringPrev=NULL, bStringNext=NULL, timeIndex=NULL, sizeFac=0.0001, NAFac=1, bStrings=NULL)
}




\arguments{
    \item{CNOlist}{
        a CNOlist structure, as created by \code{makeCNOlist}.
    }
    \item{model}{
        a model structure, as created by code{readSIF}, normally pre-processed 
        but that is not a requirement of this function.
    }
    \item{simList}{
        a simList as created by \code{prep4sim}, that has also already been cut to 
        contain only the reactions to be evaluated. If not provided, it is
        recomputed automatically.
    }
    \item{indexList}{
        a list of indexes of the species stimulated/inhibited/measured in the 
        model, as created by \code{indexFinder}. If not provided ,it is
        recomputed automatically.
    }
    \item{simResPrev}{
        Results of Previous simulation at time TN-1 step.
    }
    \item{bStringPrev}{
        the best bitString at time TN-1
    }
    \item{bStringNext}{
        the bitString to use to compute the score at time TN
    }
    \item{timeIndex}{
    Future feature will allows timeIndex to provide the exact list of indices to cut and plot. For now,
it is based on the bitStrings provided. 
    }
    \item{sizeFac}{
    the scaling factor for the size term in the objective function, default to 0.0001
    }
    \item{NAFac}{
    the scaling factor for the NA term in the objective function, default to 1    
    }
    \item{bStrings}{
    list of optimised bitstrings found at the previous time points
    }
}



\value{
    \item{score}{See \code{gaBinaryTN} for details}
}

\author{
    T. Cokelaer, S.Schrier
}

\examples{
    data(CNOlistToy2,package="CellNOptR")
    data(ToyModel2,package="CellNOptR")
    model <- preprocessing(CNOlistToy2, ToyModel2)
    bStringT1 = c(0,0,1,1,1,1,1,1,0,0,1,1,1,1,1,1)
    simT1<-simulateTN(CNOlist=CNOlistToy2, model=model, bStrings=list(bStringT1))

    score1 = computeScoreTN(CNOlistToy2, model, bStrings=list(bStringT1,c(1,0,1,0)))
}




#
#  This file is part of the CNO software
#
#  Copyright (c) 2011-2012 - EMBL - European Bioinformatics Institute
#
#  File author(s): CNO developers (cno-dev@ebi.ac.uk)
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  CNO website: http://www.cellnopt.org
#
##############################################################################
# $Id$
simulateTN<-function(CNOlist, model, bStrings){
#simResT1,bStringT2,bStringTimes,simList,indexList,timeIndex){

    if (is.list(bStrings)==FALSE){
        stop("CellNOpt Error: 3d argument called bStrings must be a list of vectors. Each vector representing a bit string")
    }

    if (length(bStrings) == 1){
        simPrev = internal_simulateT1(CNOlist, model, bStrings[[1]])
    }
    else if (length(bStrings) > 1){
        # T1 first
        simPrev = internal_simulateT1(CNOlist, model, bStrings[[1]])
        # Then, loop over T2 to TN
        for (i in 2:length(bStrings)){
            res = buildBitString(bStrings[1:i])
            simPrev = internal_simulateTN(CNOlist, model, simPrev, res$bs, res$bsTimes, i+1)
        }
    }

    return(simPrev)
}

internal_simulateTN <- function(CNOlist, model, simPrev, bStringPrev, bStringTimes, timeIndex){

  simList = prep4sim(model)
  indexList = indexFinder(CNOlist, model)
  modelCut<-cutModel(model, bStringPrev)


  simListCut<-simList
  simListCut$finalCube<-simListCut$finalCube[as.logical(bStringPrev),]
  simListCut$ixNeg<-simListCut$ixNeg[as.logical(bStringPrev),]
  simListCut$ignoreCube<-simListCut$ignoreCube[as.logical(bStringPrev),]
  simListCut$maxIx<-simListCut$maxIx[as.logical(bStringPrev)]

  if(is.null(dim(simListCut$finalCube))){
    simListCut$finalCube<-matrix(simListCut$finalCube,ncol=1)
    simListCut$ixNeg<-matrix(simListCut$ixNeg,ncol=1)
    simListCut$ignoreCube<-matrix(simListCut$ignoreCube,ncol=1)
  }
  modelCut$times <- bStringTimes[which(bStringTimes != 0)]


  # simulate
  SimT2 <- simulatorTN(
    simResultsPrev=simPrev,
    CNOlist=CNOlist,
    model=modelCut,
    simList=simListCut,
    indexList=indexList,
    timeIndex=timeIndex)

  return(SimT2)
}




internal_simulateT1<-function(CNOlist, model, bStringT1, simList=NULL, indexList=NULL){
    # simList and indexList are computed inside this function. 
    # However, for back-compatibility, we keep the arguments so that if
    # provided, we can still use them.

    # cut the model
    modelCut <- cutModel(model, bStringT1)
    if (is.null(simList)==TRUE){
        simList = prep4sim(model)
    }
    if (is.null(indexList)==TRUE){
        indexList = indexFinder(CNOlist, model)
    }

    # cut the model
    newSimList = cutSimList(simList, bStringT1)

    # compute the results
    simRes <- simulatorT1(CNOlist=CNOlist, model=modelCut, simList=newSimList, 
        indexList=indexList)

    return(simRes)
}

#
#  This file is part of the CNO software
#
#  Copyright (c) 2011-2012 - EBI
#
#  File author(s): CNO developers (cno-dev@ebi.ac.uk)
#
#  Distributed under the GPLv2 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-2.0.html
#
#  CNO website: http://www.ebi.ac.uk/saezrodriguez/software.html
#
##############################################################################
# $Id$

# Function that computes the score of a specific bitstring
# Although it is very similar to computeScoreT1, there are enough differences to
# have a different function.
computeScoreTN<-function(CNOlist, model, simList=NULL, indexList=NULL,
    simResPrev=NULL, bStringPrev=NULL, bStringNext=NULL, timeIndex=NULL,
    sizeFac=0.0001, NAFac=1, bStrings=NULL){

    # by default same behaviour as computeScoreT2
    # timeIndex=3 stands for T2 by default.
    #timeIndex = timeIndex # i.e., "tN"
    if (is.null(timeIndex)==TRUE){
        timeIndex = length(bStrings) + 1
    }

#    print(timeIndex)

    if (is.null(simList)==TRUE){
        simList = prep4sim(model)
    }
    if (is.null(indexList)==TRUE){
        indexList = indexFinder(CNOlist, model)
    }

    # if no previous results provided, we must recompute simulated results at
    # each time. This is slower but provides an easy user interface. Just
    # cnolist, model and list of bitStrings is required.
    if (is.null(simResPrev) == TRUE){
        res = buildBitString(bStrings)
        bitString = res$bs
        bStringTimes = res$bsTimes

        modelCut = cutModel(model, bitString)
        modelCut$times <- bStringTimes[which(bStringTimes != 0)]
        simListCut <- cutSimList(simList, bitString)

        if (is.null(bStrings)==FALSE){
            simResults = simulateTN(CNOlist, model, bStrings)
        }
        else{
            stop("CellNOpt erro:: you must provide either bStrings or (simResPrev, bStringPrev, bStringNext) arguments)")
        }

    }
    else{


        if (is.null(bStrings)==TRUE){
            # if no bStrings, we have prev and next bistring that must be
            # provided and it should correspond to T2 analysis
            # TODO: sanity checl that timeIndex is correct
            bStrings = list(bStringPrev, bStringNext)
        } else {
            # need to check that length of bStrings is in agreement with
            # timeIndex.
        }

        res = buildBitString(bStrings)
        bitString = res$bs
        bStringTimes = res$bsTimes

        modelCut = cutModel(model, bitString)
        modelCut$times <- bStringTimes[which(bStringTimes != 0)]
        simListCut <- cutSimList(simList, bitString)
        # Compute the simulated results
        simResults <- simulatorTN(
            simResultsPrev=simResPrev,
            CNOlist=CNOlist,
            model=modelCut,
            simList=simListCut,
            indexList=indexList,
            timeIndex=timeIndex)
    }



    #Compute the score
    Score <- getFit(
        simResults=simResults,
        CNOlist=CNOlist,
        model=modelCut,
        indexList=indexList,
        timePoint=timeIndex,
        sizeFac=sizeFac,
        NAFac=NAFac,
        nInTot=length(which(model$interMat == -1)),
		simResultsT0=NA)

    if ((class(CNOlist)=="CNOlist")==FALSE){
          CNOlist = CellNOptR::CNOlist(CNOlist)
    }

  nDataP <- sum(!is.na(CNOlist@signals[[2]]))
  Score <- Score/nDataP


  return(Score)
}

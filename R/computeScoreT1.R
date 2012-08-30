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

#Function that computes the score of a specific bitstring
# todo: this is similar to wha is done in gaBinaryT1. need to do the same for T2
computeScoreT1<-function(CNOlist, model, bString, simList=NULL, indexList=NULL, 
    sizeFac=0.0001, NAFac=1){
    # simList and indexList are computed inside this function. 
    # However, for back-compatibility, we keep the arguments so that if
    # provided, we can still use them.
    if (is.null(simList)==TRUE){
        simList = prep4sim(model)
    }
    if (is.null(indexList)==TRUE){
        indexList = indexFinder(CNOlist, model)
    }


	timeIndex = 2 # i.e., "t1"

    modelCut = cutModel(model, bString)

    simListCut <- cutSimList(simList, bString)

    # Compute the simulated results
    simResults<-simulatorT1(
        CNOlist=CNOlist,
        model=modelCut,
        simList=simListCut,
        indexList=indexList)
    # We may want to use the T0 information.
    simResultsT0<-simulatorT0(
        CNOlist=CNOlist,
        model=modelCut,
        simList=simListCut,
        indexList=indexList)

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
        simResultsT0=simResultsT0)


  if ((class(CNOlist)=="CNOlist")==FALSE){
       CNOlist = CellNOptR::CNOlist(CNOlist)
  }
  nDataP <- sum(!is.na(CNOlist@signals[[2]]))
  Score <- Score/nDataP


  return(Score)
}

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

simulatorTN <-function(simResultsPrev, CNOlist, model, simList, indexList, timeIndex=3){
    #timeIndex=3 correspond to T2

    if ((class(CNOlist)=="CNOlist")==FALSE){
        CNOlist = CellNOptR::CNOlist(CNOlist)
    }
    # check the structures
	if(is.null(CNOlist@stimuli) || is.null(CNOlist@inhibitors)) {
		stop("This function needs 'valueStimuli' and 'valueInhibitors' in CNOlist")
	}
	
	if(is.null(model$reacID) || is.null(model$namesSpecies)) {
		stop("This function needs 'reacID' and 'namesSpecies' in model")
	}

    # variables
    nStimuli <- as.integer(length(indexList$stimulated))
    nInhibitors <- as.integer(length(indexList$inhibited))
    nCond <- as.integer(dim(CNOlist@stimuli)[1])
    nReacs <- as.integer(length(model$reacID))
    nSpecies <- as.integer(length(model$namesSpecies))
    nMaxInputs <- as.integer(dim(simList$finalCube)[2])
    timeIndex <- as.integer(t(timeIndex))
    nTimes <- as.integer(length(model$times))

    # model (times vector and interMat)
    times <- as.integer(t(model$times))
    interMat <- as.integer(t(model$interMat))

    # simList
    finalCube = as.integer(as.vector(t(simList$finalCube))-1)
    ixNeg = as.integer(as.vector(t(simList$ixNeg)))
    ignoreCube = as.integer(as.vector(t(simList$ignoreCube)))
    maxIx = as.integer(simList$maxIx-1)

    # index
    indexSignals <- as.integer(as.vector(indexList$signals)-1)
    indexStimuli <- as.integer(as.vector(indexList$stimulated)-1)
    indexInhibitors <- as.integer(as.vector(indexList$inhibited)-1)
    nSignals <- length(indexSignals)

    # cnolist
    valueInhibitors <- as.integer(CNOlist@inhibitors)
    valueStimuli <- as.integer(CNOlist@stimuli)

    # simResults

    valueSimResults = as.integer(t(simResultsPrev))


    # Calling the C simulator !!
    res = .Call("simulatorTN",
        # variables
        nStimuli,
        nInhibitors,
        nCond,
        nReacs,
        nSpecies,
        nSignals,
        nMaxInputs,
        nTimes,

        # model related related
        timeIndex,
        times,
        interMat,

        # simresults
        valueSimResults, 
        # simList
        finalCube,
        ixNeg,
        ignoreCube,
        maxIx,
        # index
        indexSignals,
        indexStimuli,
        indexInhibitors,
        # cnolist
        valueInhibitors,
        valueStimuli
    )

    return(res)
}

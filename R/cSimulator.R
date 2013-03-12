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

cSimulator <- function(CNOlist, model, simList, indexList, mode=1) {

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
	
	# simList
	# used to be 
    # >>> finalCube = as.integer(as.vector(t(simList$finalCube))-1)
    # but as.vector(t is slow and can be replaced by just as.integer albeit
    # appropriate C modifications
	finalCube = as.integer(simList$finalCube-1)
	ixNeg = as.integer(simList$ixNeg)
	ignoreCube = as.integer(simList$ignoreCube)
	maxIx = as.integer(simList$maxIx-1)
	
	# index
	indexSignals <- as.integer(indexList$signals-1)
	indexStimuli <- as.integer(indexList$stimulated-1)
	indexInhibitors <- as.integer(indexList$inhibited-1)
    nSignals <- length(indexSignals)

	# cnolist
	valueInhibitors <- as.integer(CNOlist@inhibitors)
	valueStimuli <- as.integer(CNOlist@stimuli)

	res = .Call("simulatorT1",
		# variables	
		nStimuli,
		nInhibitors,
		nCond,
		nReacs,
		nSpecies,
        nSignals,
		nMaxInputs,
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
		valueStimuli,
        as.integer(mode)
	)
# should not be cut because it is used in simulateTN as an input 
#    res = res[,indexList$signals]
	return(res)
}

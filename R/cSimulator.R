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

cSimulator <- function(CNOlist, model, simList, indexList) {
	
	# check the structures
	if(is.null(CNOlist$valueStimuli) || is.null(CNOlist$valueInhibitors)) {
		stop("This function needs 'valueStimuli' and 'valueInhibitors' in CNOlist")
	}
	
	if(is.null(model$reacID) || is.null(model$namesSpecies)) {
		stop("This function needs 'reacID' and 'namesSpecies' in model")
	}
	
	# variables
	nStimuli <- as.integer(length(indexList$stimulated))
	nInhibitors <- as.integer(length(indexList$inhibited))
	nCond <- as.integer(dim(CNOlist$valueStimuli)[1])
	nReacs <- as.integer(length(model$reacID))
	nSpecies <- as.integer(length(model$namesSpecies))
	nMaxInputs <- as.integer(dim(simList$finalCube)[2])
	
	# simList
	finalCube = as.integer(as.vector(t(simList$finalCube))-1)
	ixNeg = as.integer(as.vector(t(simList$ixNeg)))
	ignoreCube = as.integer(as.vector(t(simList$ignoreCube)))
	maxIx = as.integer(simList$maxIx-1)
	
	# index
	indexSignals <- as.integer(as.vector(indexList$signals)-1)
	indexStimuli <- as.integer(as.vector(indexList$stimulated)-1)
	indexInhibitors <- as.integer(as.vector(indexList$inhibited)-1)

	# cnolist
	valueInhibitors <- as.integer(t(CNOlist$valueInhibitors))
	valueStimuli <- as.integer(t(CNOlist$valueStimuli))

	res = .Call("simulatorT1",
		# variables	
		nStimuli,
		nInhibitors,
		nCond,
		nReacs,
		nSpecies,
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
		valueStimuli
	)
	
	return(res)
}

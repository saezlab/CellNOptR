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
# $Id: simulatorT1.R 3155 2013-01-09 15:24:58Z cokelaer $

simulatorT1 <- function(CNOlist,model,simList,indexList, mode=1) {

    
	#simRes = rSimulatorT1(CNOlist, model, simList, indexList)
	simRes = cSimulator(CNOlist, model, simList, indexList, mode=mode)
	return(simRes)
}

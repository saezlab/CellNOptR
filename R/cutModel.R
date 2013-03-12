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

cutModel <- function(model, bString){
    bs = as.logical(bString)
    newmodel <- list()
    newmodel$interMat <- model$interMat[, bs]
    newmodel$notMat <- model$notMat[, bs]
    newmodel$reacID <- model$reacID[bs]
    newmodel$namesSpecies <- model$namesSpecies

    # could also add the times used in times > T1 if times
    # newmodel$times <- bStringTimes[which(bStringTimes != 0)]

    return(newmodel)
}



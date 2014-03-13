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
# $Id: checkSignals.R 4106 2013-11-06 15:53:49Z cokelaer $
checkSignals<-function(CNOlist, model ){

    if ((class(CNOlist)=="CNOlist")==FALSE){
        CNOlist = CellNOptR::CNOlist(CNOlist)
    }

    # check that model is a Model list
    if(!is.list(model)){
        stop("This function expects as input a Model as output by readSIF")
    }
    if( !all( c("reacID", "namesSpecies","interMat","notMat") %in% names(model)))  {
        stop("This function expects as input a Model as output by readSIF")
    }

    # check that all of the signals in colnames(CNOlist@signals) match to one species in model$namesSpecies
    signalsMatch<-match(colnames(CNOlist@signals[[1]]),model$namesSpecies,nomatch=0)
    if(any(signalsMatch == 0)){
        stop(paste(
            "The following signals from your CNOlist do not match any of the species in your model, and should be removed:",
            toString(colnames(CNOlist@signals[[1]])[which(signalsMatch == 0)])
            ))
    }

    # same for stimuli 
    signalsMatch<-match(colnames(CNOlist@stimuli),model$namesSpecies,nomatch=0)
    if(any(signalsMatch == 0)){
        stop(paste(
            "The following stimuli from your CNOlist do not match any of the species in your model, and should be removed:",
            toString(colnames(CNOlist@stimuli)[which(signalsMatch == 0)])
            ))
    }

    # same for inhibitors
    signalsMatch<-match(colnames(CNOlist@inhibitors),model$namesSpecies,nomatch=0)
    if(any(signalsMatch == 0)){
        stop(paste(
            "The following inhibitors from your CNOlist do not match any of the species in your model, and should be removed:",
            toString(colnames(CNOlist@inhibitors)[which(signalsMatch == 0)])
            ))
    }

}


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
# $Id: checkSignals.R 1657 2012-07-05 12:26:23Z cokelaer $
checkSignals<-function(CNOlist, model ){

    # check that CNOlist is a CNOlist
    if(!is.list(CNOlist)){
        stop("This function expects as input a CNOlist as output by makeCNOlist or normaliseCNOlist")
    }
    if(all(names(CNOlist) !=
            c("namesCues",
            "namesStimuli",
            "namesInhibitors",
            "namesSignals",
            "timeSignals",
            "valueCues",
            "valueInhibitors",
            "valueStimuli",
            "valueSignals"))){
        stop("This function expects as input a CNOlist as output by makeCNOlist")
        }

    # check that model is a Model list
    if(!is.list(model)){
        stop("This function expects as input a Model as output by readSIF")
    }
    if(all(names(model) != c("reacID", "namesSpecies","interMat","notMat"))){
        stop("This function expects as input a Model as output by readSIF")
    }

    # check that all of the signals in CNOlist$namesSignals match to one species in model$namesSpecies
    signalsMatch<-match(CNOlist$namesSignals,model$namesSpecies,nomatch=0)
    if(any(signalsMatch == 0)){
        warning(paste(
            "The following signals from your CNOlist do not match any of the species in your model, and should be removed:",
            toString(CNOlist$namesSignals[which(signalsMatch == 0)])
            ))
    }

    # same for stimuli 
    signalsMatch<-match(CNOlist$namesStimuli,model$namesSpecies,nomatch=0)
    if(any(signalsMatch == 0)){
        warning(paste(
            "The following stimuli from your CNOlist do not match any of the species in your model, and should be removed:",
            toString(CNOlist$namesStimuli[which(signalsMatch == 0)])
            ))
    }

    # same for inhibitors
    signalsMatch<-match(CNOlist$namesInhibitors,model$namesSpecies,nomatch=0)
    if(any(signalsMatch == 0)){
        warning(paste(
            "The following inhibitors from your CNOlist do not match any of the species in your model, and should be removed:",
            toString(CNOlist$namesInhibitors[which(signalsMatch == 0)])
            ))
    }

}


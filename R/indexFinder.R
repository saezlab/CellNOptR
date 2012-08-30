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
# $Id: indexFinder.R 1586 2012-06-26 14:59:24Z cokelaer $
indexFinder<-function(CNOlist, model, verbose=FALSE){

#check that CNOlist is a CNOlist
    if(!is.list(CNOlist)){
        stop("This function expects as input a CNOlist as output by makeCNOlist or normaliseCNOlist")
        }
    if(all(names(CNOlist) != c(
        "namesCues",
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

    #check that Model is a model list

    if(!is.list(model)) stop("This function expects as input a model as output by readSIF")

    if(length(model) == 4) {

        if(all(names(model) != c("reacID", "namesSpecies","interMat","notMat"))){
            stop("This function expects as input a model as output by readSIF")
            }

        }

    if(length(model) == 5) {

        if(all(names(model) != c("reacID", "namesSpecies","interMat","notMat","speciesCompressed"))){
            stop("This function expects as input a Model as output by readSIF")
            }

        }

    #Find the indexes of the signals
    signals<-match(CNOlist$namesSignals,model$namesSpecies)

    #Find the indexes of the stimulated species
    stimulated<-match(CNOlist$namesStimuli,model$namesSpecies)

    #Find the indexes of the inhibited species
    inhibited<-match(CNOlist$namesInhibitors,model$namesSpecies)

    #Print summaries
    if(verbose){
        print(paste(
            "The following species are measured:",
            toString(model$namesSpecies[signals])))
        print(paste(
            "The following species are stimulated:",
            toString(model$namesSpecies[stimulated])))
        print(paste(
            "The following species are inhibited:",
            toString(model$namesSpecies[inhibited])))
        }

    #Return a list of indexes
    indexes<-list(signals=signals,stimulated=stimulated,inhibited=inhibited)

}


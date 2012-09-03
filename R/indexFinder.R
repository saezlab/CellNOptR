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
# $Id: indexFinder.R 2293 2012-09-03 12:30:47Z cokelaer $
indexFinder<-function(CNOlist, model, verbose=FALSE){

    if ((class(CNOlist)=="CNOlist")==FALSE){
        CNOlist = CellNOptR::CNOlist(CNOlist)
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
    signals<-match(colnames(CNOlist@signals[[1]]),model$namesSpecies)

    #Find the indexes of the stimulated species
    stimulated<-match(colnames(CNOlist@stimuli),model$namesSpecies)

    #Find the indexes of the inhibited species
    inhibited<-match(colnames(CNOlist@inhibitors),model$namesSpecies)

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
    return(list(signals=signals,stimulated=stimulated,inhibited=inhibited))

}


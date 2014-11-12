#
#  This file is part of the CNO software
#
#  Copyright (c) 2011-2013 - EBI
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


# Remove the signals or cues in cnolist that are not found in model
cutCNOlist <- function(cnolist, model=NULL, cutTimeIndices=list(), verbose=FALSE){

    if (is.null(model)==TRUE && length(cutTimeIndices)==0){
        stop("Neither model nor time indices were provided. You must provide a 
model (to remove species in your cnolist that are not in the model) and/or a list 
of time indices to remove data at different time  points.")
    }

    if (is.null(model)==FALSE){
        cutCNOlist = .cutCNOlistModel(cnolist, model, verbose)
        if (length(cutTimeIndices)>0){
            cutCNOlist = .cutCNOlistTimeIndices(cutCNOlist, cutTimeIndices)
        }
    } else{ 
        if (length(cutTimeIndices)>0){
            cutCNOlist = .cutCNOlistTimeIndices(cnolist, cutTimeIndices)
        }
    }

    return(cutCNOlist)
}


.cutCNOlistModel <- function(cnolist, model, verbose=FALSE){

   cnolist_signals = colnames(cnolist@signals[[1]])
   cnolist_cues = colnames(cnolist@cues)
   species = model$namesSpecies
   times = cnolist@timepoints

   # not efficient to copy but midas files are small and this function should be
   # called only once.
   cutCNOlist = cnolist

   # Search for signals to remove
   signals2remove = c()
   for (signal in cnolist_signals){
       # if cnolist specy not in the model, we will need to remove it
       if (signal %in% species == FALSE){
           signals2remove = c(signals2remove, signal)
       }
   }

   # Search for cues to remove
   cues2remove = c()
   for (cue in cnolist_cues){
       if (cue %in% species == FALSE){
           cues2remove = c(cues2remove, cue)
       }
   }

   if (verbose == TRUE){
       cat("The following signals and cues were not found in the model\n")
       cat("They will be removed from the cnolist\n")
       cat("signals:", signals2remove, "\n")
       cat("cues:",cues2remove, "\n")
   }

   # let us work with the indices instead of the names
   indices_signals = c()
   for (signal in signals2remove){
        # note that in the grep function with use the ^ and $ sign to ensure
        # that a pattern signal matches exactly a signal (e.g., STAT does not
        # match STAT1)
        indices_signals = c(indices_signals,
            grep(paste("^", signal,  "$", sep=""), cnolist_signals))
   }

   indices_cues = c()
   for (cue in cues2remove){
        indices_cues = c(indices_cues,
            grep(paste("^", cue, "$", sep=""), cnolist_cues))
   }

   # now perform the cut on the cues
   if (length(indices_cues)>0){
       cutCNOlist@cues = t(as.matrix(cnolist@cues[,-indices_cues]))
   }

   # and the signal matrices
   if (length(indices_signals)>0){
       for (time in 1:length(cnolist@signals)){
           cutCNOlist@signals[[time]] = as.matrix(cnolist@signals[[time]][,-indices_signals])
       }
    }

   return(cutCNOlist)
}
        


.cutCNOlistTimeIndices <- function(cnolist, cutTimeIndices){

    cutCNOlist = cnolist
    indices = unlist(cutTimeIndices)
    if (max(indices)>length(cnolist@signals)){
        warning("one or more indices provided are larger than number of time
        points. These indices are going to be ignored.")
    }
   cutCNOlist@signals = cnolist@signals[-indices]
   cutCNOlist@timepoints = cnolist@timepoints[-indices]

   return(cutCNOlist)
}

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


# Remove the signals, sstimulis, inhibitors or cues in cnolist that are not found in model
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
  
  cnolist_stimuli = colnames(cnolist@stimuli)
  cnolist_signals = colnames(cnolist@signals[[1]])
  cnolist_cues = colnames(cnolist@cues)
  cnolist_inhibitors = colnames(cnolist@inhibitors)
  species = model$namesSpecies
  times = cnolist@timepoints
  
  # not efficient to copy but midas files are small and this function should be
  # called only once.
  cutCNOlist = cnolist
  
  # Search for stimuli to remove
  stimuli2remove = c()
  for(stimuli in cnolist_stimuli){
    if(stimuli %in% species == FALSE){
      stimuli2remove = c(stimuli2remove, stimuli)
    }
  }
  
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
  
  # Search for inhibitors to remove
  inhibitors2remove = c()
  for (inhibitor in cnolist_inhibitors){
    if (inhibitor %in% species == FALSE){
      inhibitors2remove = c(inhibitors2remove, inhibitor)
    }
  }
  
  if (verbose == TRUE){
    cat("The following signals, cues and inhibitors were not found in the model\n")
    cat("They will be removed from the cnolist\n")
    cat("signals:", signals2remove, "\n")
    cat("cues:",cues2remove, "\n")
    cat("inhibitors:", inhibitors2remove, "\n")
  }
  
  # let us work with the indices instead of the names
  indeces_stimuli = c()
  for(stimuli in stimuli2remove){
    indeces_stimuli = c(indeces_stimuli,
                        grep(paste("^", stimuli, "$", sep = ""), cnolist_stimuli))
  }
  
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
  
  indices_inhibitors = c()
  for (inhibitor in inhibitors2remove){
    indices_inhibitors = c(indices_inhibitors,
                           grep(paste("^", inhibitor, "$", sep=""), cnolist_inhibitors))
  }
  
  # now perform cut on stimuli
  if(length(indeces_stimuli)>0){
    #cutCNOlist@stimuli = as.matrix(cnolist@stimuli[, -indeces_stimuli])
  	cutCNOlist@stimuli = cnolist@stimuli[, -indeces_stimuli,drop=FALSE]   # rewrote to keep the names of the matrixes
  }
  
  # now perform the cut on the cues
  if (length(indices_cues)>0){
    #cutCNOlist@cues = as.matrix(cnolist@cues[,-indices_cues])
  	cutCNOlist@cues = cnolist@cues[,-indices_cues,drop=FALSE]
  }
  
  # and the signal matrices
  if (length(indices_signals)>0){
    for (time in 1:length(cnolist@signals)){
      #cutCNOlist@signals[[time]] = as.matrix(cnolist@signals[[time]][,-indices_signals])
      #cutCNOlist@variances[[time]] = as.matrix(cnolist@variances[[time]][,-indices_signals])
      cutCNOlist@signals[[time]] = cnolist@signals[[time]][,-indices_signals,drop=FALSE]
      cutCNOlist@variances[[time]] = cnolist@variances[[time]][,-indices_signals,drop=FALSE]
    }
  }
  
  # and the inhibitor matrices
  if (length(indices_inhibitors)>0){
    #for (time in 1:length(cnolist@inhibitors)){
    #cutCNOlist@inhibitors[[time]] = as.matrix(cnolist@inhibitors[[time]][,-indices_inhibitors])
    #}
    #cutCNOlist@inhibitors = as.matrix(cnolist@inhibitors[,-indices_inhibitors])
    cutCNOlist@inhibitors = cnolist@inhibitors[,-indices_inhibitors,drop=FALSE]
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

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
#  CNO website: http://www.ebi.ac.uk/saezrodriguez/cno
#
##############################################################################
# $Id$
prep4sim<-function(model){
  
  #Compute the max number of inputs observed in the model for a single reaction
  maxInput<-colSums(model$interMat)
  maxInput<-abs(min(maxInput))+1
  
  #Make the empty matrices
  finalCube<-matrix(NA, nrow=length(model$reacID),ncol=maxInput)
  ixNeg<-matrix(FALSE, nrow=length(model$reacID),ncol=maxInput)
  ignoreCube<-matrix(TRUE, nrow=length(model$reacID),ncol=maxInput)
  maxIx<-rep(NA,length(model$reacID))
  
  #Fill the matrices finalCube, ignoreCube and ixNeg, and maxIx
  
  for(r in 1:length(model$reacID)){
    
    input<-which(model$interMat[,r] == -1)
      
    if (length(input)>0){
      finalCube[r,1:length(input)]<-input
      
      if(length(input) < maxInput) finalCube[r,(length(input)+1):maxInput]<-1
      
      neg<-model$notMat[input,r]
      ixNeg[r,1:length(input)]<-(neg == 1)
      ignoreCube[r,1:length(input)]<-FALSE
      maxIx[r]<-which(model$interMat[,r] == 1)
    }
    
  }
  
  rownames(finalCube)<-model$reacID
  rownames(ixNeg)<-model$reacID
  rownames(ignoreCube)<-model$reacID
  names(maxIx)<-model$reacID
  
  modelname<-deparse(substitute(model))
  
  # check for self-loops here
  # rescan reactions for self-loops
  for(r in 1:length(model$reacID)){
    # pick out species
    lhs = strsplit(model$reacID[r], "=")[[1]][1]
    rhs = strsplit(model$reacID[r], "=")[[1]][2]
    lhs = unlist(strsplit(lhs, "\\+"))
    # notmat is okay so no need to record if negative or not
    lhs = gsub("!","",lhs)
    if(any(lhs==rhs)) { # this is a self loop with/without AND gate
        input = sapply(lhs, function(x) which(model$namesSpecies == x))
        
        finalCube[r,1:length(input)]<-input
        if(length(input) < maxInput) finalCube[r,(length(input)+1):maxInput]<-1
        
        neg<-model$notMat[input,r]
        ixNeg[r,1:length(input)]<-(neg == 1)
        ignoreCube[r,1:length(input)]<-FALSE
        maxIx[r]<-which(model$interMat[,r] == 1)
        
      }
    }
  
  Fields4Sim<-list(
    finalCube=finalCube,
    ixNeg=ixNeg,
    ignoreCube=ignoreCube,
    maxIx=maxIx,
    modelname=modelname, maxInput=maxInput)
  
  return(Fields4Sim)
}

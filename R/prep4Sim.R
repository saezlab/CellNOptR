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
#  CNO website: http://www.ebi.ac.uk/saezrodriguez/cno
#
##############################################################################
# $Id: prep4Sim.R 595 2012-02-22 17:21:47Z cokelaer $
prep4Sim<-function(Model){

#Compute the max number of inputs observed in the model for a single reaction
	maxInput<-colSums(Model$interMat)
	maxInput<-abs(min(maxInput))+1
	
#Make the empty matrices
	finalCube<-matrix(NA, nrow=length(Model$reacID),ncol=maxInput)
	ixNeg<-matrix(FALSE, nrow=length(Model$reacID),ncol=maxInput)
	ignoreCube<-matrix(TRUE, nrow=length(Model$reacID),ncol=maxInput)
	maxIx<-rep(NA,length(Model$reacID))
	
#Fill the matrices finalCube, ignoreCube and ixNeg, and maxIx

	for(r in 1:length(Model$reacID)){
	
		input<-which(Model$interMat[,r] == -1)
		finalCube[r,1:length(input)]<-input
		
		if(length(input) < maxInput) finalCube[r,(length(input)+1):maxInput]<-1	
		
		neg<-Model$notMat[input,r]
		ixNeg[r,1:length(input)]<-(neg == 1)
		ignoreCube[r,1:length(input)]<-FALSE
		maxIx[r]<-which(Model$interMat[,r] == 1)
		
		}
		
	rownames(finalCube)<-Model$reacID
	rownames(ixNeg)<-Model$reacID
	rownames(ignoreCube)<-Model$reacID
	names(maxIx)<-Model$reacID
	
	modelname<-deparse(substitute(Model)) 
	
	Fields4Sim<-list(
		finalCube=finalCube,
		ixNeg=ixNeg,
		ignoreCube=ignoreCube,
		maxIx=maxIx,
		modelname=modelname)
	
	return(Fields4Sim)
	}


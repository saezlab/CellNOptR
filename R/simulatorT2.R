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
# $Id$
simulatorT2 <-function(
	SimResultst1,
	CNOlist,
	Model,
	SimList,
	indexList){
	
	nSp<-dim(Model$interMat)[1]
	nReacs<-dim(Model$interMat)[2]
	nCond<-dim(CNOlist$valueStimuli)[1]
	maxIpg<-dim(SimList$finalCube)[2]
	
	if(is.null(dim(Model$interMat))){ 
		nSp<-length(Model$interMat)
		nReacs<-1
		maxIpg<-length(SimList$finalCube)
	}
	
	#This holds, for each sp, how many reactions have that sp as output
	endIx<-rep(NA,nSp)
	for(i in 1:nSp){
		endIx[i]<-length(which(SimList$maxIx == i))
		}
	maxgpo<-max(endIx)	
	
#This value is used to test the stop condition for difference between 2 iterations	
	testVal<-1E-3
	
#Create an initial values matrix	
	initValues<-SimResultst1
	
#Initialise main loop
	newInput<-initValues
	termCheck1<-TRUE
	termCheck2<-TRUE
	count<-1
	
#First iteration
	outputPrev<-newInput
	tempStore<-apply(SimList$finalCube,2,function(x){return(outputPrev[,x])})
	
		filltempCube<-function(x){
			cMatrix<-matrix(data=x,nrow=dim(SimList$ixNeg)[1],ncol=nCond)
			cVector<-apply(cMatrix,1,function(x){return(x)})
			return(cVector)
			}
			
	tempIxNeg<-apply(SimList$ixNeg,2,filltempCube)
	tempIgnore<-apply(SimList$ignoreCube,2,filltempCube)
	tempStore[tempIgnore]<-NA
	tempStore[tempIxNeg]<-1-tempStore[tempIxNeg]
	
	minNA<-function(x){
		if(all(is.na(x))){
			return(NA)
			}else{
				return(min(x,na.rm=TRUE))
				}
		}
		
	outputCube<-apply(tempStore,1,minNA)
	outputCube<-matrix(outputCube, nrow=nCond,ncol=nReacs)
	
		#rewrite anything that comes to a node that also receives a t2 branch,ie set it to the same 
		#as our t2 reac for those reacs, so that they won't influence the OR
	reacsT2<-which(Model$times == 2)
	
	if(length(reacsT2) > 0){
	
		for(i in reacsT2){
			outNode<-which(Model$interMat[,i] > 0)
			reacs2Overwrite<-which(Model$interMat[outNode,] > 0)
			
			if(length(reacs2Overwrite) != 0) {
				for(n in 1:length(reacs2Overwrite)){
					outputCube[,reacs2Overwrite[n]]<-outputCube[,i]
					}
				}
				
			}
			
		}	
		
	for(s in 1:nSp){
	
		if(endIx[s] != 0){
			compOR<-function(x){
				if(all(is.na(x[which(SimList$maxIx == s)]))){
					res<-NA
					}else{
						res<-max(x[which(SimList$maxIx == s)],na.rm=TRUE)
						}
				return(res)
				}
			newInput[,s]<-apply(outputCube,1,compOR)
			}
			
		}
		
	for(stim in 1:length(indexList$stimulated)){
	
		stimM<-cbind(CNOlist$valueStimuli[,stim],newInput[,indexList$stimulated[stim]])
		
		maxNA<-function(x){
			return(max(x,na.rm=TRUE))
			}
			
		stimV<-apply(stimM,1,maxNA)
		newInput[,indexList$stimulated[stim]]<-stimV
		
		}
		
	valueInhibitors<-1-CNOlist$valueInhibitors
	newInput[,indexList$inhibited]<-valueInhibitors*newInput[,indexList$inhibited]
	newInput[is.na(newInput)]<-0
	outputPrev[is.na(outputPrev)]<-0
	termCheck1<-!all(abs(outputPrev-newInput)<testVal)
	termCheck2<-(count < (nSp*1.2))
	count<-count+1
	firstIter<-newInput
	
#Main loop

	while(termCheck1 && termCheck2){
	
		outputPrev<-newInput
		
#This is now a 2 columns matrix that has a column for each input (column in finalCube)
#and a set of rows for each reac (where a set contains as many rows as conditions)
#all concatenated into one long column
		tempStore<-apply(SimList$finalCube,2,function(x){return(outputPrev[,x])})
		tempIxNeg<-apply(SimList$ixNeg,2,filltempCube)
		tempIgnore<-apply(SimList$ignoreCube,2,filltempCube)
		
#Set to NA the values that are "dummies", so they won't influence the min		
		tempStore[tempIgnore]<-NA
		
#Flip the values that enter with a negative sign		
		tempStore[tempIxNeg]<-1-tempStore[tempIxNeg]
		
#Compute all the ands by taking, for each gate, the min value across the inputs of that gate			
		outputCube<-apply(tempStore,1,minNA)
		
#OutputCube is now a vector of length (nCond*nReacs) that contains the input of each reaction in
#each condition, concatenated as such allcond4reac1,allcond4reac2,etc...		
#This is transformed into a matrix with a column for each reac and a row for each cond		
		outputCube<-matrix(outputCube, nrow=nCond,ncol=nReacs)
		
#Go through each species, and if it has inputs, then take the max across those input reactions
#i.e. compute the ORs
		for(s in 1:nSp){
			if(endIx[s] != 0){
				compOR<-function(x){
					if(all(is.na(x[which(SimList$maxIx == s)]))){
						res<-NA
						}else{
							res<-max(x[which(SimList$maxIx == s)],na.rm=TRUE)
							}
					return(res)
					}
				newInput[,s]<-apply(outputCube,1,compOR)
				}
			}
			
#Reset the inhibitors and stimuli
		for(stim in 1:length(indexList$stimulated)){
			stimM<-cbind(CNOlist$valueStimuli[,stim],newInput[,indexList$stimulated[stim]])
			maxNA<-function(x){
				return(max(x,na.rm=TRUE))
				}
			stimV<-apply(stimM,1,maxNA)
			newInput[,indexList$stimulated[stim]]<-stimV
			}
		valueInhibitors<-1-CNOlist$valueInhibitors
		newInput[,indexList$inhibited]<-valueInhibitors*newInput[,indexList$inhibited]

#Set all the nodes that are targets of a t2 reaction to the state that they had at the first iteration
		if(length(reacsT2) != 0){
			t2reacs<-Model$interMat[,reacsT2]
			for(r in 1:length(reacsT2)){
				if(length(reacsT2) == 1) {
					target<-which(t2reacs > 0)
					}else{
						target<-which(t2reacs[,r] > 0)
						}
				newInput[,target]=firstIter[,target]
				}
			}
			
#Replace NAs with zeros to avoid having the NA	penalty applying to unconnected species
		newInput[is.na(newInput)]<-0
		outputPrev[is.na(outputPrev)]<-0

#Check the 2 termination conditions		
		termCheck1<-!all(abs(outputPrev-newInput)<testVal)
		termCheck2<-(count < (nSp*1.2))
		count<-count+1
		}
		
#Set the non-resolved bits to NA
	newInput[which(abs(outputPrev-newInput) > testVal)]<-NA
	return(newInput)
	
	}


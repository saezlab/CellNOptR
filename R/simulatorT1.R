simulatorT1<-function(CNOlist,Model,SimList,indexList){
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
	initValues<-matrix(data=NA,nrow=nCond,ncol=nSp)
	colnames(initValues)<-Model$namesSpecies
	
#Set the initial values of the stimuli	
	initValues[,indexList$stimulated]<-CNOlist$valueStimuli
	
#Flip the inhibitors so that 0=inhibited/1=non-inhibitted	
	valueInhibitors<-1-CNOlist$valueInhibitors
	valueInhibitors[which(valueInhibitors == 1)]<-NA
	
#Set the initial values of the inhibited species: 0 if inhibited, untouched if not inhibited	
	initValues[,indexList$inhibited]<-valueInhibitors
	
#Initialise main loop
	newInput<-initValues
	termCheck1<-TRUE
	termCheck2<-TRUE
	count<-1
	
#Main loop
	while(termCheck1 && termCheck2){
	
		outputPrev<-newInput
#This is now a 2 columns matrix that has a column for each input (column in finalCube)
#and a set of rows for each reac (where a set contains as many rows as conditions)
#all concatenated into one long column

		filltempCube<-function(x){
			cMatrix<-matrix(data=x,nrow=nReacs,ncol=nCond)
			cVector<-apply(cMatrix,1,function(x){return(x)})
			return(cVector)
			}
			
		if(nReacs > 1){
		
			tempStore<-apply(SimList$finalCube,2,function(x){return(outputPrev[,x])})	
			tempIxNeg<-apply(SimList$ixNeg,2,filltempCube)
			tempIgnore<-apply(SimList$ignoreCube,2,filltempCube)
			
			}else{
			
				tempStore<-outputPrev[,SimList$finalCube]
				tempIxNeg<-matrix(
					SimList$ixNeg,nrow=nCond,
					ncol=length(SimList$ixNeg),byrow=TRUE)
				tempIgnore<-matrix(
					SimList$ignoreCube,nrow=nCond,
					ncol=length(SimList$ignoreCube),byrow=TRUE)
				
				}
#Set to NA the values that are "dummies", so they won't influence the min		
		tempStore[tempIgnore]<-NA
		
#Flip the values that enter with a negative sign		
		tempStore[tempIxNeg]<-1-tempStore[tempIxNeg]
		
		minNA<-function(x){
			if(all(is.na(x))){
				return(NA)
				}else{
					return(min(x,na.rm=TRUE))
					}
			}
			
#Compute all the ands by taking, for each gate, the min value across the inputs of that gate			
		if(nReacs > 1){
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
				
			}else{
				outputCube<-ifelse(all(is.na(tempStore)),NA,min(tempStore,na.rm=TRUE))
				newInput[,SimList$maxIx]<-outputCube
				}
				
#Reset the inhibitors and stimuli
		for(stim in 1:length(indexList$stimulated)){
			stimM<-cbind(
				CNOlist$valueStimuli[,stim],
				newInput[,indexList$stimulated[stim]])
			maxNA<-function(x){
				return(max(x,na.rm=TRUE))
				}
				
			stimV<-apply(stimM,1,maxNA)
			newInput[,indexList$stimulated[stim]]<-stimV
			}
			
		valueInhibitors<-1-CNOlist$valueInhibitors
		newInput[,indexList$inhibited]<-valueInhibitors*newInput[,indexList$inhibited]
		
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


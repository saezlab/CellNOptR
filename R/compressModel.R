compressModel<-function(Model, indexes){

#this will hold the indices of the species that we have to look at (i.e. not signal nor cue)
	species<-seq(1:length(Model$namesSpecies))
	species<-species[-c(indexes$signals,indexes$stimulated,indexes$inhibited)]
	speciesCompressed<-0
	
	for(sp in species){
	
#Find all the reacs that come out of species sp
		outReac<-which(Model$interMat[sp,] == -1)
		
#Find all the reacs that come into species sp
		inReac<-which(Model$interMat[sp,] == 1)
		
#If the species doesn't have any output, we can't compress it, so we do sthg only if
#the species has one output minimum
		if(length(outReac) != 0){
		
#If the species has only one input		
			if(length(inReac) == 1){
			
#Find the input element(s) (there might be more than one element if the input is an & gate)
				input<-which(Model$interMat[,inReac] == -1)
				
#Find the output elements
				if(length(outReac) > 1) 
					outputs<-apply(Model$interMat[,outReac],2,function(x) which(x == 1))
				if(length(outReac) == 1) 
					outputs<-which(Model$interMat[,outReac] == 1)
				
#check that the input is not an & gate, and that none of the downstream elements is the output element				
				if(length(input) == 1 && !any(outputs == input[1])){
				
#need to create new reacs, that go from input to outputs, and to remove the old reacs

#Create the first one, that goes from input to the first output
					reac2Remove<-c(inReac,outReac)
					newReacs<-Model$interMat[,inReac]+Model$interMat[,outReac[1]]
					newNots<-Model$notMat[,inReac]
					if(Model$notMat[sp,outReac[1]] == 1){
						newNots[input]<-1
						}
					newReacIDs<-paste(
						names(newReacs)[which(newReacs == -1)],
						"=",
						names(newReacs)[which(newReacs == 1)],
						sep="")
					if(any(newNots == 1)){
						newReacIDs<-paste("!",newReacIDs,sep="")
						}
#now create the ones that go to the following outputs, if any
					if(length(outReac) > 1){
					
						for(r in 2:length(outReac)){
							newReacs<-cbind(
								newReacs,
								(Model$interMat[,inReac]+Model$interMat[,outReac[r]]))
							nNtemp<-Model$notMat[,inReac]
							if(Model$notMat[sp,outReac[r]] == 1) nNtemp[input]<-1
							newNots<-cbind(newNots,nNtemp)
							newReacIDstemp<-paste(
								rownames(newReacs)[which(newReacs[,r] == -1)],
								"=",
								rownames(newReacs)[which(newReacs[,r] == 1)],
								sep="")
							if(any(newNots[,r] == 1)){
								newReacIDstemp<-paste("!",newReacIDstemp,sep="")
								}
							newReacIDs<-c(newReacIDs,newReacIDstemp)
							}
							
						colnames(newReacs)<-newReacIDs
						colnames(newNots)<-newReacIDs
						}
						
#now remove the reacs to remove and appende the ones that have just been created
					Model$interMat<-Model$interMat[,-reac2Remove]
					Model$notMat<-Model$notMat[,-reac2Remove]
					Model$reacID<-Model$reacID[-reac2Remove]
					Model$interMat<-cbind(Model$interMat, newReacs)
					Model$notMat<-cbind(Model$notMat,newNots)
					Model$reacID<-c(Model$reacID,newReacIDs)
					
					if(length(outReac) == 1){
						colnames(Model$interMat)[dim(Model$interMat)[2]]<-newReacIDs
						colnames(Model$notMat)[dim(Model$notMat)[2]]<-newReacIDs
						}
						
#Store which species has been compressed. This will be used to remove rows of the matrices
#after all is compressed (I don't do it know because then the indexes of the species that we
#have decided to look at in the beginning won't be valid anymore)
					speciesCompressed<-c(speciesCompressed,sp)
					}
#if: check that the input is not an & gate, 
#and that none of the downstream elements is the output element	
#then the above was done -> else:
					
				}else{
				
#If the species has more than one input, but only one output			
					if(length(outReac) == 1){
					
#Find the output element
						output<-which(Model$interMat[,outReac] == 1)
						
#Find the input element(s)
						if(length(inReac) > 1){
							inputs<-apply(Model$interMat[,inReac],2,function(x) which(x == -1))
							}
							
						if(length(inReac) == 1){
							inputs<-which(Model$interMat[,inReac] == -1)
							}
						
#check that none of the downstream elements is the output element	
						if(!any(inputs == output[1])){
						
#create the first reac, that goes from the first input to the output
							reac2Remove<-c(inReac,outReac)
							newReacs<-Model$interMat[,inReac[1]]+Model$interMat[,outReac]
							newNots<-Model$notMat[,inReac[1]]
							
#if the in gate is not an &, but the out is negated:the in is inverted	(if it was neg it becomes pos, and vice versa)				
							if(Model$notMat[sp,outReac] == 1 && length(which(newReacs == -1)) == 1){
								newNots[inputs[1]]<-ifelse((newNots[inputs[1]] == 1), 0, 1)
								}

#if the in gate is an &, and the out is negated:
							if(Model$notMat[sp,outReac] == 1 && length(which(newReacs == -1)) != 1) {
								newNots[inputs[1]]<-ifelse((newNots[inputs[1]] == 1), 0, 1)
								newNots[inputs[2]]<-ifelse((newNots[inputs[2]] == 1), 0, 1)
								}

#This makes the name of the new reac in cases where the input is not an & gate							
							if(length(which(newReacs == -1)) == 1) {
								newReacIDs<-paste(
									names(newReacs)[which(newReacs == -1)],
									"=",
									names(newReacs)[which(newReacs == 1)],
									sep="")
								if(any(newNots == 1)) newReacIDs<-paste("!",newReacIDs,sep="")
								indone<-1
								}else{

#This makes the name of the new reac in cases where the input is an & gate								
									in1<-names(newReacs)[which(newReacs == -1)[1]]
									in2<-names(newReacs)[which(newReacs == -1)[2]]
									if(newNots[which(newReacs == -1)[1]] == 1) in1<-paste("!",in1,sep="")
									if(newNots[which(newReacs == -1)[2]] == 1) in2<-paste("!",in2,sep="")
									newReacIDs<-paste(
										in1,"+",in2,"=",
										names(newReacs)[which(newReacs == 1)],
										sep="")
									indone<-2
									}

#Now create the reactions that go to from other inputs	to the output

							if(length(inReac) > 1){
							
								for(r in 2:length(inReac)){
								
									newReacs<-cbind(newReacs, (Model$interMat[,inReac[r]]+Model$interMat[,outReac]))
									nNtemp<-Model$notMat[,inReac[r]]
									
									if(Model$notMat[sp,outReac] == 1 && length(which(newReacs == -1)) == 1){
										nNtemp[inputs[indone+1]]<-ifelse((nNtemp[inputs[indone+1]] == 1), 0, 1)
										}
										
									if(Model$notMat[sp,outReac] == 1 && length(which(newReacs == -1)) != 1) {
										nNtemp[inputs[indone+1]]<-ifelse((nNtemp[inputs[indone+1]] == 1), 0, 1)
										nNtemp[inputs[indone+2]]<-ifelse((nNtemp[inputs[indone+2]] == 1), 0, 1)
										}
										
									newNots<-cbind(newNots,nNtemp)	
									
									if(length(which(newReacs[,dim(newReacs)[2]] == -1)) == 1){
									
										newReacIDstemp<-paste(
											rownames(newReacs)[which(newReacs[,r] == -1)],
											"=",
											rownames(newReacs)[which(newReacs[,r] == 1)],
											sep="")
											
										if(any(newNots[,r] == 1)){
											newReacIDstemp<-paste("!",newReacIDstemp,sep="")
											}
										newReacIDs<-c(newReacIDs,newReacIDstemp)
										indone<-indone+1
										
										}else{
										
											in1<-rownames(newReacs)[which(newReacs[,r] == -1)[1]]
											in2<-rownames(newReacs)[which(newReacs[,r] == -1)[2]]
											
											if(newNots[which(newReacs == -1)[1],r] == 1) in1<-paste("!",in1,sep="")
											
											if(newNots[which(newReacs == -1)[2],r] == 1) in2<-paste("!",in2,sep="")
											
											newReacIDstemp<-paste(
												in1,"+",in2,"=",
												rownames(newReacs)[which(newReacs[,r] == 1)],
												sep="")
											newReacIDs<-c(newReacIDs,newReacIDstemp)
											indone<-indone+2
											}
											
									colnames(newReacs)<-newReacIDs
									colnames(newNots)<-newReacIDs
									
									}
								}
								
							Model$interMat<-Model$interMat[,-reac2Remove]
							Model$notMat<-Model$notMat[,-reac2Remove]
							Model$reacID<-Model$reacID[-reac2Remove]
							Model$interMat<-cbind(Model$interMat, newReacs)
							Model$notMat<-cbind(Model$notMat,newNots)
							Model$reacID<-c(Model$reacID,newReacIDs)
							
							if(length(inReac) == 1){
								colnames(Model$interMat)[dim(Model$interMat)[2]]<-newReacIDs
								colnames(Model$notMat)[dim(Model$notMat)[2]]<-newReacIDs
								}	
								
							speciesCompressed<-c(speciesCompressed,sp)	
							}
						}
					}

#If the species has more than one input and more than one output then nothing is done		
			}
		}
	ModelCompressed<-Model
	ModelCompressed$speciesCompressed<-NA

#If species have been compressed
	if(length(speciesCompressed) > 1){

#Remove the dummy 0 in speciesCompressed		
		speciesCompressed<-speciesCompressed[2:length(speciesCompressed)]
		ModelCompressed$speciesCompressed<-Model$namesSpecies[speciesCompressed]
		ModelCompressed$namesSpecies<-Model$namesSpecies[-speciesCompressed]
		ModelCompressed$notMat<-Model$notMat[-speciesCompressed,]
		ModelCompressed$interMat<-Model$interMat[-speciesCompressed,]
		}

#Check that no duplicates reactions have been created
	uniqueR<-unique(ModelCompressed$reacID)
	
	if(length(uniqueR) != length(ModelCompressed$reacID)){
		ModelCompressed$reacID<-ModelCompressed$reacID[match(uniqueR,ModelCompressed$reacID)]
		ModelCompressed$interMat<-ModelCompressed$interMat[,match(uniqueR,colnames(ModelCompressed$interMat))]
		ModelCompressed$notMat<-ModelCompressed$notMat[,match(uniqueR,colnames(ModelCompressed$notMat))]
		}
		
	return(ModelCompressed)	
	
	}


cutNONC<-function(Model,NONCindexes){

	if(length(NONCindexes) == 0){
	
		newModel=Model
		
		}else{
		
			newspecies<-Model$namesSpecies[-NONCindexes]
			newInterMat<-Model$interMat[-NONCindexes,]
			newNotMat<-Model$notMat[-NONCindexes,]
			
#this function finds whether a given vector contains at least one input and one output or not
#it returns true if the re is an in/out missing
			EmptyInOut<-function(x){
			
				input<-match(-1,x,nomatch=0)
				output<-match(1,x,nomatch=0)
				
				if((input == 0) | (output == 0)){
					return(TRUE)
					}else{
						return(FALSE)
						}
						
				}
				
			reac2remove<-apply(newInterMat,2,EmptyInOut)
			
			if(any(reac2remove)){
			
				reac2remove<-which(reac2remove)
				newInterMat<-newInterMat[,-reac2remove]
				newNotMat<-newNotMat[,-reac2remove]
				newreacID<-Model$reacID[-reac2remove]
				
				}else{
				
					newModel<-list(
						reacID=Model$reacID,
						namesSpecies=newspecies,
						interMat=newInterMat,
						notMat=newNotMat)
					
					}
			
			newModel<-list(
				reacID=newreacID,
				namesSpecies=newspecies,
				interMat=newInterMat,
				notMat=newNotMat)
			
			}
			
	return(newModel)
	}


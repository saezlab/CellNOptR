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


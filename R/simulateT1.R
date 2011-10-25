simulateT1<-function(CNOlist,Model,bStringT1,SimList,indexList){

	Modelcut<-Model
	
	Modelcut$interMat<-Modelcut$interMat[,as.logical(bStringT1)]
	Modelcut$notMat<-Modelcut$notMat[,as.logical(bStringT1)]
	Modelcut$reacID<-Modelcut$reacID[as.logical(bStringT1)]
	
	SimListcut<-SimList
	SimListcut$finalCube<-SimListcut$finalCube[as.logical(bStringT1),]
	SimListcut$ixNeg<-SimListcut$ixNeg[as.logical(bStringT1),]
	SimListcut$ignoreCube<-SimListcut$ignoreCube[as.logical(bStringT1),]
	SimListcut$maxIx<-SimListcut$maxIx[as.logical(bStringT1)]
	
	if(is.null(dim(SimListcut$finalCube))){
		SimListcut$finalCube<-matrix(SimListcut$finalCube,ncol=1)
		SimListcut$ixNeg<-matrix(SimListcut$ixNeg,ncol=1)
		SimListcut$ignoreCube<-matrix(SimListcut$ignoreCube,ncol=1)
		}
		
	SimRes<-simulatorT1(
		CNOlist=CNOlist,
		Model=Modelcut,
		SimList=SimListcut,
		indexList=indexList)
	
	return(SimRes)
	}


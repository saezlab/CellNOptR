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
# $Id: simulateT1.R 595 2012-02-22 17:21:47Z cokelaer $
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


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
model2sif<-function(model,optimRes=NA,writeSif=FALSE, filename="Model"){
	
  if (is.na(optimRes[1])){
    BStimes<-rep(1,length(model$reacID))
  }else{
    BStimes<-optimRes$bString
  }	

	findOutput<-function(x){
		sp<-which(x == 1)
		sp<-model$namesSpecies[sp]
		}
		
	reacOutput<-apply(model$interMat,2,findOutput)
	
	findInput<-function(x){
		sp<-which(x == -1)
		sp<-model$namesSpecies[sp]
		}
		
	reacInput<-apply(model$interMat,2,findInput)
		
#if the class of reacInput is not a list, then there are no AND gates
	if(class(reacInput) != "list"){
	
		isNeg<-function(x){
			isNegI<-any(x == 1)
			return(isNegI)
			}
			
		inpSign<-apply(model$notMat,2,isNeg)
		inpSign<-!inpSign
		inpSign[inpSign]<-1
		inpSign[!inpSign]<--1
    
		sifFile<-cbind(reacInput,inpSign,reacOutput)
    sifFile<-sifFile[BStimes==1,]
		colnames(sifFile)=NULL
    rownames(sifFile)=NULL
    
		}else{
		
#in this case there are AND gates and so we need to create dummy "and#" nodes
			sifFile<-matrix(ncol=3)
			nANDs<-1
			for(i in 1:length(reacOutput)){
			  if (BStimes[i]==1){

				  if(length(reacInput[[i]]) == 1){
            tmp<-matrix(0,nrow=1,ncol=3)
					  tmp[1,1]<-reacInput[[i]]
					  tmp[1,3]<-reacOutput[i]
					  tmp[1,2]<-ifelse(
						  any(model$notMat[,i] == 1),-1,1)
              sifFile<-rbind(sifFile,tmp)
          
					}else{
					
						for(inp in 1:length(reacInput[[i]])){
              tmp<-matrix(0,nrow=1,ncol=3)
							tmp[1,1]<-reacInput[[i]][inp]
							tmp[1,3]<-paste("and",nANDs,sep="")
							tmp[1,2]<-ifelse(
								model$notMat[which(reacInput[[i]][inp]==rownames(model$notMat)),i] == 1,
								-1,1)
              sifFile<-rbind(sifFile,tmp)
							}
						tmp<-matrix(0,nrow=1,ncol=3)	
						tmp[1,1]<-paste("and",nANDs,sep="")
						tmp[1,3]<-reacOutput[i]
						tmp[1,2]<-1
            sifFile<-rbind(sifFile,tmp)
						
            nANDs<-nANDs+1
            }
			    }
				}
				
      sifFile<-sifFile[2:dim(sifFile)[1],]
			
			}

#this is the edge attribute matrix that contains, for each edge, whether it is
#absent from the model (0), present at t1(1) or present at t2(2)
  if (writeSif==TRUE){
	filename<-paste(filename, ".sif", sep="")
    writeSIF(sifFile, filename=filename)
  }

  return(sifFile)
}

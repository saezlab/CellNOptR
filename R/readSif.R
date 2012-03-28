# -*- python -*-
#
#  This file is part of the CNO software
#
#  Copyright (c) 2011-2012 - EBI
#
#  File author(s): Thomas Cokelaer <cokelaer@ebi.ac.uk>
#
#  Distributed under the GPLv2 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-2.0.html
#
#  CNO website: http://www.ebi.ac.uk/saezrodriguez/software.html
#
##############################################################################
# $Id: readSif.R 591 2012-02-22 17:17:14Z cokelaer $
readSif<-function(sifFile){

#Read the sif file
	sif<-read.table(sifFile)
	sif<-as.matrix(sif)
	
#Create the vector of names of species	
	namesSpecies<-unique(c(as.character(sif[,1]), as.character(sif[,3])))
	
#Create the vector of names or reactions	

	createReacID<-function(x){
		neg<-ifelse(test=(x[2] == "-1"), TRUE, FALSE)
		ID<-ifelse(
			neg, 
			paste("!",x[1],"=",x[3],sep=""), 
			paste(x[1],"=",x[3],sep=""))
		return(ID)
		}
		
	reacID<-apply(sif,1,createReacID)
	
#Create empty interMat and notMat	
	interMat<-matrix(data=0,nrow=length(namesSpecies), ncol=dim(sif)[1])
	rownames(interMat)<-namesSpecies
	colnames(interMat)<-reacID
	notMat<-matrix(data=0,nrow=length(namesSpecies), ncol=dim(sif)[1])
	rownames(notMat)<-namesSpecies
	colnames(notMat)<-reacID
	
#Fill in interMat and notMat	
	for(r in 1:length(reacID)){
	
		source<-match(sif[r,1],namesSpecies)
		target<-match(sif[r,3],namesSpecies)
		neg<-ifelse(test=(sif[r,2] == "-1"), TRUE, FALSE)
		interMat[source,r]<-(-1)
		interMat[target,r]<-1
		if(neg) notMat[source,r]<-1
		
		}
		
#Take care of the "and" nodes

#1.detect them
	AndNodes<-grep(
		pattern="(and\\d+$)",namesSpecies,perl=TRUE,ignore.case=FALSE)
	AndNodesV<-grep(
		pattern="(and\\d+$)",namesSpecies,perl=TRUE,ignore.case=FALSE,value=TRUE)
	if(length(AndNodes) != 0){
	
#2.go through them	
		for(i in 1:length(AndNodes)){
		
#3.store the name of the node, and find the associated reactions		
			node<-AndNodesV[i]
			AndsReacs<-c(
				grep(pattern=paste(node,"=",sep=""),reacID,ignore.case=FALSE),
				grep(pattern=paste(node,"$",sep=""),perl=TRUE,reacID,ignore.case=FALSE))

#4.create the new column of interMat and notMat			
			newReac<-rowSums(interMat[,AndsReacs])
			newReacNot<-rowSums(notMat[,AndsReacs])
			
#5.create the new reacID

			#First case: I have an "and" with two inputs only
			
			if(length(AndsReacs) == 3){
			
				output<-grep(pattern=paste(node,"=",sep=""),reacID,
					ignore.case=FALSE,value=TRUE)
				inputs<-grep(pattern=paste(node,"$",sep=""),perl=TRUE,
					reacID,ignore.case=FALSE,value=TRUE)
				output<-sub(pattern=node,x=output,replacement="",ignore.case=FALSE)
				inputs<-sub(pattern=paste("=",node,sep=""),x=inputs,
					replacement="+",ignore.case=FALSE)
				inputs<-paste(inputs,collapse='')
				newReacID<-paste(inputs,output,sep="")
				newReacID<-sub(pattern=paste("\\W{1}",output,sep=""),
					x=newReacID,replacement=output,ignore.case=FALSE,perl=TRUE)
				
				}else{
				
			#Second case: I have an 'and' with more than 2 inputs
					output<-grep(pattern=paste(node,"=",sep=""),reacID,
						ignore.case=FALSE,value=TRUE)
						
					if(length(output) == 1){
					
						output<-sub(pattern=node,x=output,replacement="",ignore.case=FALSE)
						inputs<-grep(pattern=paste(node,"$",sep=""),perl=TRUE,
							reacID,ignore.case=FALSE,value=TRUE)
						inputs<-sub(pattern=paste("=",node,sep=""),x=inputs,
							replacement="+",ignore.case=FALSE)
						inputs<-paste(inputs,collapse='')
						newReacID<-paste(inputs,output,sep="")
						newReacID<-sub(pattern=paste("\\W{1}",output,sep=""),
							x=newReacID,replacement=output,ignore.case=FALSE,perl=TRUE)
						
						}else{
						
							output<-sub(pattern=paste(node,"=",sep=""),x=output,
								replacement="",ignore.case=FALSE)
							inputs<-grep(pattern=paste(node,"$",sep=""),perl=TRUE,
								reacID,ignore.case=FALSE,value=TRUE)
							inputs<-sub(pattern=paste("=",node,sep=""),x=inputs,
								replacement="+",ignore.case=FALSE)
							inputs<-paste(inputs,collapse='')
							inputs<-sub(pattern="\\W{1}$",x=inputs,replacement="",
								ignore.case=FALSE,perl=TRUE)
							newReacID<-paste(inputs,output,sep="=")
							
							}
					}
					
#5b.if there are more than one outputs, then the interMat and notMat created above are not correct
			output<-grep(pattern=paste(node,"=",sep=""),reacID,
				ignore.case=FALSE,value=FALSE)
			inputs<-grep(pattern=paste(node,"$",sep=""),perl=TRUE,
				reacID,ignore.case=FALSE,value=FALSE)
				
			if(length(output) != 1){
			
				newReac<-rowSums(interMat[,c(inputs,output[1])])
				newReacNot<-rowSums(notMat[,c(inputs,output[1])])
				
				for(n in 2:length(output)){
					newReac<-cbind(newReac,rowSums(interMat[,c(inputs,output[n])]))
					newReacNot<-cbind(newReacNot,rowSums(notMat[,c(inputs,output[n])]))
					}
					
				}
				
#6. remove the old reactions and replace them by the new ones
			namesSpecies<-namesSpecies[-which(namesSpecies == AndNodesV[i])]
			reacID<-reacID[-AndsReacs]
			interMat<-interMat[,-AndsReacs]
			notMat<-notMat[,-AndsReacs]
			interMat<-cbind(interMat,newReac)
			notMat<-cbind(notMat,newReacNot)
			reacID<-c(reacID,newReacID)
			interMat<-interMat[-which(rownames(interMat) == AndNodesV[i]),]
			notMat<-notMat[-which(rownames(notMat) == AndNodesV[i]),]
			
			if(length(output) == 1){
			
				colnames(interMat)[dim(interMat)[2]]<-newReacID
				colnames(notMat)[dim(notMat)[2]]<-newReacID
				
				}else{
				
					colnames(interMat)[(dim(interMat)[2]-length(newReacID)+1):dim(interMat)[2]]<-newReacID
					colnames(notMat)[(dim(notMat)[2]-length(newReacID)+1):dim(notMat)[2]]<-newReacID
					
					}
			}
		}
		
	Model<-list(
		reacID=reacID,
		namesSpecies=namesSpecies,
		interMat=interMat,
		notMat=notMat)	
	return(Model)
	}


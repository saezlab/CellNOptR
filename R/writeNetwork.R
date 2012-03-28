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
# $Id: writeNetwork.R 592 2012-02-22 17:18:16Z cokelaer $
##the global function is still called writeNetwork and still needs the same arguments 
#(ModelComprExpanded, optimResT1, optimResT2, ModelOriginal, CNOlist), but it is divided into: 
#	-a function writeSNetworkW that does the actual writing to files, 
#	-a function that gets the info for the sif and dot files
writeNetwork<-function(
	ModelOriginal,
	ModelComprExpanded,
	optimResT1,
	optimResT2=NA,
	CNOlist, 
    tag=NULL,verbose=FALSE){

	nwInfo<-getNetworkInfo(
		ModelOriginal=ModelOriginal,
		ModelComprExpanded=ModelComprExpanded,
		optimResT1=optimResT1,
		optimResT2=optimResT2,
		CNOlist=CNOlist,verbose=verbose)
		
	writeNetworkW(
		dN=nwInfo$dN,
		dM=nwInfo$dM,
		ModelOriginal=ModelOriginal,
		sifFile=nwInfo$sifFile,
		EAweights=nwInfo$EAweights,
		nodesAttr=nwInfo$nodesAttr, 
        tag=tag)
	
	}
	
	
#########################################################################
##############these are the functions used above
#########################################################################
#this is the bit that writes
writeNetworkW<-function(dN,dM,ModelOriginal,sifFile,EAweights,nodesAttr, tag=NULL){

 create_filename<-function(x, tag=NULL){
        if (is.null(tag)){
            return(x)
        }
        else{
            return(paste(c(tag, "_", x), collapse=""))
        }
    }


	writeDot(
		dotNodes=dN,
		dotMatrix=dM,
		Model=ModelOriginal,
		filename=create_filename("PKN.dot", tag=tag))
		
	write.table(
		sifFile[,1:3],
		file=create_filename("PKN.sif", tag=tag),
		row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
		
	write.table(
		EAweights,
		file=create_filename("TimesPKN.EA", tag=tag),
		row.names=FALSE,col.names="Times",quote=FALSE,sep="\t")	
		
	write.table(
		nodesAttr,
		file=create_filename("nodesPKN.NA", tag=tag),
		row.names=FALSE,col.names="NodesType",quote=FALSE,sep="\t")	
		
	}	
############
#this is the bit that produces the info that is needed above
getNetworkInfo<-function(
	ModelOriginal,
	ModelComprExpanded,
	optimResT1,
	optimResT2,
	CNOlist,
	verbose){

	#These are used to create the string that holds the information about whether an edge is present
#or not (BStimes)
	bString1<-optimResT1$bString
	
	if(is.na(optimResT2[1])){
		bString2<-optimResT1$bString[which(optimResT1$bString == 0)]
		}else{
			bString2<-optimResT2$bString
			}
			
#BStimes is a string containing 0,1 and 2 depending on whether the interaction is 
#absent, present at t1 or present at t2
	BStimes<-bString1
	BStimes[which(BStimes == 0)]<-bString2*2
#!!! BStimes relates to strings in the compressed and expanded model
	reacs<-ModelComprExpanded$reacID
	reacsOriginal<-ModelOriginal$reacID
	
#check if there are ANDs that were split, and if yes, set the score of the original AND 
#to the max score of its resulting branches (i.e.the max of the branches, if one branch was present
#at t1 the AND gets 1, same if one was present at t2, if one of the branches was present at T1
#and another one at t2, then the AND gets a score of 3)

	if(names(ModelComprExpanded$SplitANDs) != "initialReac"){
	
		for(i in 1:length(ModelComprExpanded$SplitANDs)){
		
			map<-match(ModelComprExpanded$SplitANDs[[i]],reacs)
			mapScore<-BStimes[map]
		#this forces the score to be 3 if a branch was present at t1 and another at t2	
			if(any(mapScore == 1) && any(mapScore == 3)) mapScore[1]<-3
			mapScore<-max(mapScore)
			splitAnd<-names(ModelComprExpanded$SplitANDs)[i]
			splitAnd<-which(reacs == splitAnd)
			BStimes[splitAnd]<-mapScore
			
			}
			
		}
		
#check if there were ands that were created, and if yes, assign to each branch of the original OR
#the max of the score of the and or of that particular branch

	if(length(ModelComprExpanded$newANDs) != 0){
	
		for(i in 1:length(ModelComprExpanded$newANDs)){
			branches<-ModelComprExpanded$newANDs[i]
			branchesPos<-match(unlist(branches),reacs)
			BStimes[unlist(branches)]<-pmax(
				BStimes[match(unlist(branches),reacs)],
				BStimes[match(names(ModelComprExpanded$newANDs)[i],reacs)])
			}
			
		}		
		
#The new BStimes vector now contain weights that have been adapted to propagate weights coming 
#from new ANDs or split ANDs
#Now create an adjacency list for the compressed/expanded model: col1=in node, col2=out node, col3=weight

findOutput<-function(x){
		sp<-which(x == 1)
		sp<-ModelComprExpanded$namesSpecies[sp]
		}
		
findInput<-function(x){
		sp<-which(x == -1)
		sp<-ModelComprExpanded$namesSpecies[sp]
		}
		
n<-1
adj<-matrix(NA,nrow=(length(reacs)*3),ncol=3)

for(i in 1:length(reacs)){

	input<-findInput(ModelComprExpanded$interMat[,i])
	output<-findOutput(ModelComprExpanded$interMat[,i])
	weight<-BStimes[i]
	
	for(j in 1:length(input)){
		adj[n,1]<-input[j]
		adj[n,2]<-output
		adj[n,3]<-weight
		n<-n+1
		}
		
	}
	
adj<-adj[-which(is.na(adj[,1])),]

#remove duplicates that arise from new ANDs and split ANDs
adj<-cbind(
	adj[order(adj[,3],decreasing=TRUE),1],
	adj[order(adj[,3],decreasing=TRUE),2],
	adj[order(adj[,3],decreasing=TRUE),3])
adj<-unique(adj)

for(i in 1:dim(unique(adj[,1:2]))[1]){

	findMatch<-apply(adj,1,function(x) all(x[1:2] == adj[i,1:2]))
	
	if(sum(findMatch) > 1){
		findMatch[which(findMatch == TRUE)[1]]<-FALSE
		adj<-adj[-which(findMatch==TRUE),]
		}
		
	}
	
#Now do the same for the original model: col1=in node, col2=out node, col3=reaction #

findOutput<-function(x){
		sp<-which(x == 1)
		sp<-ModelOriginal$namesSpecies[sp]
		}
		
findInput<-function(x){
		sp<-which(x == -1)
		sp<-ModelOriginal$namesSpecies[sp]
		}
		
n<-1
adjOrig<-matrix(NA,nrow=(length(reacsOriginal)*3),ncol=3)

for(i in 1:length(reacsOriginal)){

	input<-findInput(ModelOriginal$interMat[,i])
	output<-findOutput(ModelOriginal$interMat[,i])
	
	for(j in 1:length(input)){
		adjOrig[n,1]<-input[j]
		adjOrig[n,2]<-output
		adjOrig[n,3]<-i
		n<-n+1
		}
		
	}
	
adjOrig<-adjOrig[-which(is.na(adjOrig[,1])),]

#Go through the adjacency list of the compressed model and find which path is the shortest 
#between those species
#this will hold the max weight of anything that has been mapped to this row of adjOrig so far
OrigMap<-rep(0,dim(adjOrig)[1])

for(i in 1:dim(adj)[1]){
#If the row of the compressed adjacency list matches a row of the original adjacency list,
#we just copy the weights
	rowAdj<-adj[i,1:2]
	matchR<-apply(adjOrig[,1:2],1,function(x) all(x == rowAdj))
	
	if(sum(matchR) > 0){
	
		OrigMap[which(matchR)]<-max(OrigMap[which(matchR)],adj[i,3])
		
		}else{
		
#Otherwise we need to find the shortest path between the nodes in that row of adj, within the rows
#of adjOrig, and map the weights of the row of adj to each of the row of adjOrig that are involved 
#in that path
#This only works for paths of length 2 at most: improve this!!!!!
			ins<-which(adjOrig[,1] == adj[i,1])
			outN<-adjOrig[ins,2]
			target<-FALSE
			n<-1
			
			while(n <= length(outN) && !target){
				outs<-which(adjOrig[,1] == outN[n])
				outN2<-adjOrig[outs,2]
				
				if(any(outN2 == adj[i,2])){
				
					rOut<-which(outN2 == adj[i,2])
					rOut<-outs[rOut]
					rIn<-ins[n]
					n<-n+1
					target<-TRUE
					
					}else{
					
						n<-n+1
						#this sets those variables to NA if the path couldn't be mapped
						rIn<-NA
						rOut<-NA
						
						}
						
				}	
			#if the path could be mapped, we record the mapping	
			if(!is.na(rIn) && !is.na(rOut)){
			
				OrigMap[rIn]<-max(OrigMap[rIn],adj[i,3])	
				OrigMap[rOut]<-max(OrigMap[rOut],adj[i,3])
			#if the path couldn't be mapped (path too long), we print a warning, nad OrigMap 
			#stays 0
				}else{
					
					if(verbose == TRUE){
						print("Please be aware that when mapping the scaffold network back to the PKN, compressed paths of length > 2 are ignored.")
						}
					
					}
			
			}
	}
	BStimesOrig<-cbind(adjOrig[,3],OrigMap)
	
	for(i in 1:dim(BStimesOrig)[1]){
		BStimesOrig[i,2]<-max(BStimesOrig[i,2],BStimesOrig[which(BStimesOrig[,2] == BStimesOrig[i,1]),2])
		}
		
	BStimesOrig<-unique(BStimesOrig)	
	BStimesOrig<-BStimesOrig[,2]
	
#These mini functions are used to find the outputs and inputs of a reaction

	findOutput<-function(x){
		sp<-which(x == 1)
		sp<-ModelOriginal$namesSpecies[sp]
		}
		
	reacOutput<-apply(ModelOriginal$interMat,2,findOutput)
	
	findInput<-function(x){
		sp<-which(x == -1)
		sp<-ModelOriginal$namesSpecies[sp]
		}
		
	reacInput<-apply(ModelOriginal$interMat,2,findInput)
	
#This mini function is used to create a reaction label as used in a cystoscape edge attribute file

	createReac<-function(x){
		r<-paste(x[1]," (",x[2],") ",x[3],sep="")
		return(r)
		}
		
#if the class of reacInput is not a list, then there are no AND gates

	if(class(reacInput) != "list"){
	
		isNeg<-function(x){
			isNegI<-any(x == 1)
			return(isNegI)
			}
			
		inpSign<-apply(ModelOriginal$notMat,2,isNeg)
		inpSign<-!inpSign
		inpSign[inpSign]<-1
		inpSign[!inpSign]<--1
		sifFile<-cbind(reacInput,inpSign,reacOutput)
		EApresent<-apply(sifFile,1,createReac)
		EAweights<-cbind(EApresent, BStimesOrig)
		
		}else{
		
#in this case there are AND gates and so we need to create dummy "and#" nodes
			sifFile<-matrix(0,nrow=200,ncol=4)
			nR<-1
			nANDs<-1
			
			for(i in 1:length(reacOutput)){
			
				if(length(reacInput[[i]]) == 1){
				
					sifFile[nR,1]<-reacInput[[i]]
					sifFile[nR,3]<-reacOutput[i]
					sifFile[nR,2]<-ifelse(any(ModelOriginal$notMat[,i] == 1),-1,1) 
					sifFile[nR,4]<-BStimesOrig[i]
					nR<-nR+1
					
					}else{
					
						for(inp in 1:length(reacInput[[i]])){
						
							sifFile[nR,1]<-reacInput[[i]][inp]
							sifFile[nR,3]<-paste("and",nANDs,sep="")
							sifFile[nR,2]<-ifelse(ModelOriginal$notMat[inp,i] == 1,-1,1)
							sifFile[nR,4]<-BStimesOrig[i]
							nR<-nR+1
							
							}
							
						sifFile[nR,1]<-paste("and",nANDs,sep="")
						sifFile[nR,3]<-reacOutput[i]
						sifFile[nR,2]<-1
						sifFile[nR,4]<-BStimesOrig[i]
						nANDs<-nANDs+1	
						nR<-nR+1
						
						}
				}
				
			sifFile<-sifFile[1:(nR-1),]
			EApresent<-apply(sifFile[,1:3],1,createReac)
			EAweights<-cbind(EApresent,sifFile[,4])
			
			}
	#this bit creates a matrix of reaction that will be used for the dot file
	dM<-cbind(sifFile,BStimesOrig)
	
	#this mini function makes edge attributes in the format e1 (sign) e2 = attr		
	
	makeEA<-function(x){
		ea<-paste(x[1],"=",x[2])
		return(ea)
		}
		
	EAweights<-apply(EAweights,1,makeEA)	
	
#write the nodes attributes
	nodesCompr<-ModelComprExpanded$speciesCompressed
	indexes<-indexFinder(CNOlist,ModelOriginal)
	nodesSig<-ModelOriginal$namesSpecies[indexes$signals]
	nodesInh<-ModelOriginal$namesSpecies[indexes$inhibited]
	nodesStim<-ModelOriginal$namesSpecies[indexes$stimulated]
	nodesNCNO<-findNONC(ModelOriginal,indexes)
	nodesAttrNames<-nodesSig
	nodesAttr<-rep("signal",length(nodesSig))
	
	if(length(nodesInh) != 0){
		nodesAttrNames<-c(nodesAttrNames,nodesInh)
		nodesAttr<-c(nodesAttr,rep("inhibited",length(nodesInh)))
		}
		
	if(length(nodesStim) != 0){
		nodesAttrNames<-c(nodesAttrNames,nodesStim)
		nodesAttr<-c(nodesAttr,rep("stimulated",length(nodesStim)))
		}	
		
	if(length(nodesNCNO) != 0){
		nodesAttrNames<-c(nodesAttrNames,nodesNCNO)
		nodesAttr<-c(nodesAttr,rep("ncno",length(nodesNCNO)))
		}		
		
	if(length(nodesCompr) != 0){
		nodesAttrNames<-c(nodesAttrNames,nodesCompr)
		nodesAttr<-c(nodesAttr,rep("compressed",length(nodesCompr)))
		}	
		
	nodesAttr<-cbind(nodesAttrNames,nodesAttr)	
	
#this bit creates a matrix of nodes attributes that will be used for the dot file
	dN<-nodesAttr
	
#this is the node attributes in the format as will be used by cytoscape	
	nodesAttr<-apply(nodesAttr,1,makeEA)	
	
#nodes that are neither inh/stim/sig nor NONC will not be in thi snode attribute so I'll add them
	other<-setdiff(ModelOriginal$namesSpecies,dN[,1])
	other<-paste(other," = NA",sep="")
	nodesAttr<-c(nodesAttr, other)
	
#this is the edge attribute matrix that contains, for each edge, whether it is
#absent from the model (0), present at t1(1) or present at t2(2)
	return(list(
		dN=dN,
		dM=dM,
		sifFile=sifFile,
		nodesAttr=nodesAttr,
		EAweights=EAweights))
	}

findNONC<-function(Model,indexes,verbose=FALSE){

#use the floyd warshall algorithm implemented in the package RBGL
#eg floyd.warshall.all.pairs.sp(coex) that returns a distance matrix
#for that I need to create a graphNEL object

	ModeltoGraphNEL<-function(Model){
	
		vertices<-Model$namesSpecies
		edgeList<-vector("list", length=length(vertices))
		names(edgeList)<-vertices	
	#edgeList is a list with an element edges and an element weights, both of the same size
	#where each element	in edgeL refers to one node of the graph and contains a vector made of
	#the indexes of the nodes to which the node in question has an edge
	
		for(sp in 1:length(vertices)){
			reacs<-which(Model$interMat[sp,] == -1)
			
			if(length(reacs) == 0){
				targets<-sp
				
				}else{	
					targets<-which(Model$interMat[,reacs[1]] == 1)
					
					if(length(reacs) > 1){
					
						for(r in 2:length(reacs)){
							targets<-c(targets, which(Model$interMat[,reacs[r]] == 1))
							}
							
						}
					}
	#If the node doesn't have a target, I set it a dummy interaction with itself because
	#otherwise the graphNEL object can't be created. This won't change whether there exists 
	#a path to it from a signal/cue or not
			edgeList[[sp]] <- list(edges=targets, weights=rep(1,length(targets)))
			}
			
		graph<-new("graphNEL", nodes=vertices, edgeL=edgeList,edgemode="directed")
		return(graph)
		}
		
#Now I use the function above to create a GraphNEL object		
	graphModel<-ModeltoGraphNEL(Model)	
	
#now I can run the floyd warshall algo
	distMatrix<-floyd.warshall.all.pairs.sp(graphModel)
	
#Find the species that are non-observable, i.e. there are no path from that species to any
#of the signals
#And the species that are not controllable, ie there are no paths from any of the cues
#to those species
	ncno<-0
	
	for(sp in 1:length(Model$namesSpecies)){
#if all the distances from this species to signals are infinite and the species is not a signal	
#then this species is not observable
		if(
			all(is.infinite(distMatrix[sp,indexes$signals]))  && 
			!any(indexes$signals == sp)){
			
			ncno<-c(ncno,sp)
			
			}
			
		if(
			all(is.infinite(distMatrix[c(indexes$stimulated,indexes$inhibited),sp])) && 
			!any(c(indexes$stimulated,indexes$inhibited) == sp)){
			
			ncno<-c(ncno,sp)
			
			}
		}
	ncno<-ncno[-1]
	ncno<-unique(ncno)
	
#check that none of the cues or signals are in there
	cs<-intersect(unlist(indexes),ncno)
	
	if(length(cs != 0)){
		ncno<-ncno[-match(cs,ncno)]
		}
		
	if(verbose){
		print(paste(
			"The following species are not observable and/or not controllable:",
			toString(Model$namesSpecies[ncno])))
		}
	return(ncno)
	}


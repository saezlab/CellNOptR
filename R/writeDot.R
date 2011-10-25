writeDot<-function(dotNodes,dotMatrix,Model,fileName){

#Write the dot file for the PKN:
#Use the internal variables created by writeNetwork:
	#dotNodes is a matrix with 2 columns: the first has the node names, 
		#and the second the attributes (signal, stimulated, inhibited, compressed, ncno)
		#! a node can appear twice in this matrix if it belongs to more of one of the above categories
		#! a node could also not appear here if it is is none of these categories
	#dotMatrix is a matrix with 4 or 5 columns, and a row for each reaction:
		#the first column holds the name of the input node, 
		#the second column holds the sign of the reaction (-1 if negative, 1 if positive)
		#the third column holds the name of the ouptut node
		#an optional 5th column holds the weights of the edges
		#the fourth column holds the time stamp (0,1,2)
		
#1. Create a vector of nodes and their attribute (class: stimulated, inhibited , measured,etc.)
	nodes<-unique(c(dotMatrix[,1],dotMatrix[,3]))	
	nClass<-rep(NA,length(nodes))
	
#this creates a vector of nodes category (stimulated, compressed, etc) in the same order as
#the nodes vector.  If a node is not in any of these categories, the entry will stay NA.

	for(i in 1:length(nodes)){
	#match only returns the first match, which implies the following order of priority for the
	#colors: signals prime over inhib which primes over stim which primes over ncno which primes 
	#over compressed
	#if I later want to put different colors on the nodes with more than one attribute,
	#I could still use, instead of match, which(dotNodes[,1] == nodes[1]), which would give me 
	#a vector when the node has more than one attribute, I could then look into the identity 
	#of the entries of this vector and define a color depending on this
		nClass[i]<-dotNodes[match(nodes[i],dotNodes[,1]),2]
		}
		
	AndNodes<-grep(pattern="(and\\d+$)",nodes,perl=TRUE,ignore.case=FALSE)	
	
	if(length(AndNodes) != 0){
		nClass[AndNodes]<-"dummy"
		}
		
#2.Find the ranks of nodes
#source nodes = those that do not get any input
#sink nodes = those that do not have any output

	findSources<-function(x){
		notSource<-any(x == 1)
		return(notSource)
		}
		
	notSources<-apply(Model$interMat,1,findSources)
	sources<-rownames(Model$interMat[!notSources,])
	
	findSinks<-function(x){
		notSink<-any(x == -1)
		return(notSink)
		}
		
	notSinks<-apply(Model$interMat,1,findSinks)
	sinks<-rownames(Model$interMat[!notSinks,])
	
#now find the distances from nodes to source/sink
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
		
	graphModel<-ModeltoGraphNEL(Model)	
	distMatrix<-floyd.warshall.all.pairs.sp(graphModel)
	rankNodes<-rep(Inf, length(nodes))
	
	for (i in 1:length(nodes)){
		classNode<-ifelse(is.na(nClass[i]),"na",nClass[i])
	
	if(classNode != "dummy"){
			rankNodes[i]<-min(distMatrix[sources,nodes[i]])
		
		if(any(sinks == nodes[i])){
				rankNodes[i]<-"sink"
				
				}else{
				
					if(any(sources == nodes[i])){
						rankNodes[i]<-"source"
						}
					}
			}		
		}
#Correct the ranks for 'and' nodes (I'll just put them at the same level as their inputs)
	nClass[is.na(nClass)]<-"na"
	
	if(any(nClass == "dummy")){
		dummyNodes<-nodes[which(nClass == "dummy")]
		dummyRank<-rep(1,length(dummyNodes))
		
		for(i in length(dummyNodes)){
			inN<-dotMatrix[which(dotMatrix[,3] == dummyNodes[i]),1]
			ranksIn<-rankNodes[match(inN,nodes)]
			dummyRank[i]<-min(ranksIn)
		
		}
		rankNodes[which(nClass == "dummy")]<-dummyRank

		}

	nClass[which(nClass == "na")]<-NA
	
#3.Write the header and nodes part of the dot file
#Write the header
	cat('digraph G{\nsize="8.5,11";\n{rank=source;',file=fileName,append=TRUE,sep="")
	cat(which(rankNodes == "source"),file=fileName,append=TRUE,sep=";")
	cat(';}\n{rank=same;',file=fileName,append=TRUE,sep="")
	cat(which(rankNodes == "1"),file=fileName,append=TRUE,sep=";")
	
	for(i in 2:max(rankNodes[-which(rankNodes=="source"|rankNodes=="sink")])){
		cat(';}\n{rank=same;',file=fileName,append=TRUE,sep="")
		cat(which(rankNodes == i),file=fileName,append=TRUE,sep=";")
		}
		
	cat(';}\n',file=fileName,append=TRUE,sep="")
	cat('{rank=sink;',file=fileName,append=TRUE,sep="")
	cat(which(rankNodes == "sink"),file=fileName,append=TRUE,sep=";")
	cat(';}\n',file=fileName,append=TRUE,sep="")
	
#write the nodes
	for(i in 1:length(nodes)){
		cat(i,file=fileName,append=TRUE,sep="")
		
		if(is.na(nClass[i])){
		
			cat(' [color="black" shape="ellipse" style="solid" label="',
				file=fileName,append=TRUE,sep="")
			cat(nodes[i],file=fileName,append=TRUE,sep="")
			cat('" fontname=Helvetica fontsize=22.0 ];\n',
				file=fileName,append=TRUE,sep="")
		
			}else{
			
				if(nClass[i] == "signal"){
					cat(' [color="lightblue" shape="ellipse" style="filled" label="',
						file=fileName,append=TRUE,sep="")
					}
					
				if(nClass[i] == "inhibited"){
					cat(' [color="orangered" shape="ellipse" style="filled" label="',
						file=fileName,append=TRUE,sep="")
					}
					
				if(nClass[i] == "stimulated"){
					cat(' [color="olivedrab3" shape="ellipse" style="filled" label="',
						file=fileName,append=TRUE,sep="")
					}
					
				if(nClass[i] == "ncno"){
					cat(' [color="lightblue" shape="ellipse" style="filled" label="',
						file=fileName,append=TRUE,sep="")
					}
					
				if(nClass[i] == "compressed"){
					cat(' [color="black" shape="ellipse" style="dashed" label="',
						file=fileName,append=TRUE,sep="")
					}
					
				if(nClass[i] == "dummy"){
					cat(' [color="grey90" shape="circle" style="filled" label="" fontname=Helvetica fontsize=22.0 fixedsize=true width="0.050000" height="0.050000" ];\n',file=fileName,append=TRUE,sep="")
					}
					
				if(nClass[i] != "dummy"){
					cat(nodes[i],file=fileName,append=TRUE,sep="")
					cat('" fontname=Helvetica fontsize=22.0 ];\n',
						file=fileName,append=TRUE,sep="")
					}	
					
				}
		}
#4.Write the reaction part of the dot file
#in dotMatrix, the column 2 determines the arrowhead (neg/pos), the 4 determines the color 
#grey90 if edge not in optimal model, forestgreen if at t1, and blue if at t1 and/or t2

	if(dim(dotMatrix)[2] < 5){
	
		weightsE<-rep(4,dim(dotMatrix)[1])
		
		}else{
		
			weightsE<-1+(as.numeric(dotMatrix[,5])*3)
			
			}
			
#here I need to insert the bit to make the weights from 1 to 4
	for (i in 1:dim(dotMatrix)[1]){
	
#a. get the reaction input->output
		cat(
			paste(match(dotMatrix[i,1],nodes),"->",match(dotMatrix[i,3],nodes),sep=" "),
			file=fileName,append=TRUE,sep="")
		
#b.get the color
		if(dotMatrix[i,4] == 0){
			cat('[ color="grey90" label="" weight="1.000000" penwidth="',
				file=fileName,append=TRUE,sep="")
			}
			
		if(dotMatrix[i,4] == 1){
			cat('[ color="forestgreen" label="" weight="1.000000" penwidth="',
				file=fileName,append=TRUE,sep="")
			}	
			
		if(dotMatrix[i,4] == 2){
			cat('[ color="blue" label="" weight="1.000000" penwidth="',
				file=fileName,append=TRUE,sep="")
			}	
			
#c. get the penwidth from  the weights	
		cat(weightsE[i],file=fileName,append=TRUE,sep="")
		
#get the arrowhead
		if(dotMatrix[i,2] == 1){
			cat('" arrowhead="normal" style="solid"];\n',file=fileName,append=TRUE,sep="")
			}
			
		if(dotMatrix[i,2] == -1){
			cat('" arrowhead="tee" style="solid"];\n',file=fileName,append=TRUE,sep="")
			}	
		}
	cat("}",file=fileName,append=TRUE,sep="")	
}
	
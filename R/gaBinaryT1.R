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
# $Id: gaBinaryT1.R 804 2012-03-22 16:56:26Z cokelaer $
gaBinaryT1<-function(
	CNOlist,
	Model,
	SimList,
	indexList,
	sizeFac=0.0001,
	NAFac=1,
	initBstring,
	PopSize=50,
	Pmutation=0.5,
	MaxTime=60,
	maxGens=500,
	StallGenMax=100,
	SelPress=1.2,
	elitism=5, 
	RelTol=0.1, 
	verbose=TRUE){
	
#initialise
	bLength<-length(initBstring)
	Pop<-rbind(
		initBstring,
		round(matrix(runif(bLength*(PopSize-1)), nrow=(PopSize-1),ncol=bLength)))
	bestbit<-Pop[1,]
	bestobj<-Inf
	stop<-FALSE
	obj<-rep(0,PopSize)
	g<-0
	stallGen<-0
	res<-rbind(
		c(g,bestobj,toString(bestbit),stallGen,Inf,Inf,toString(bestbit),0),
		c(g,bestobj,toString(bestbit),stallGen,Inf,Inf,toString(bestbit),0))
	colnames(res)<-c(
		"Generation",
		"Best_score",
		"Best_bitString",
		"Stall_Generation",
		"Avg_Score_Gen",
		"Best_score_Gen",
		"Best_bit_Gen",
		"Iter_time")
	PopTol<-rep(NA,bLength)
	PopTolScores<-NA
	
#Function that produces the score for a specific bitstring
	getObj<-function(x){
	
		bitString<-x
		
	#cut the model according to bitstring	
		ModelCut<-Model
		ModelCut$interMat<-ModelCut$interMat[,as.logical(bitString)]
		ModelCut$notMat<-ModelCut$notMat[,as.logical(bitString)]
		ModelCut$reacID<-ModelCut$reacID[as.logical(bitString)]

		SimListCut<-cutSimList(SimList, bitString)
		
	#compute the simulated results	
		SimResults<-simulatorT1(
			CNOlist=CNOlist,
			Model=ModelCut,
			SimList=SimListCut,
			indexList=indexList)

    # We may want to to use the T0 information.
		SimResultsT0<-simulatorT0(
			CNOlist=CNOlist,
			Model=ModelCut,
			SimList=SimListCut,
			indexList=indexList)
		
	#Compute the score	
		Score<-getFit(
			SimResults=SimResults,
			SimResultsT0=SimResultsT0,
			CNOlist=CNOlist,
			Model=ModelCut,
			indexList=indexList,
			timePoint="t1",
			sizeFac=sizeFac,
			NAFac=NAFac,
			nInTot=length(which(Model$interMat == -1)))
		nDataP<-sum(!is.na(CNOlist$valueSignals[[2]]))
		Score<-Score/nDataP

		return(Score)
		}
		
#Loop
	t0<-Sys.time()
	t<-t0
	
	while(!stop){
	
		#compute the scores
		scores<-apply(Pop,1,getObj)
		
		#Fitness assignment: ranking, linear
		rankP<-order(scores,decreasing=TRUE)
		Pop<-Pop[rankP,]
		scores<-scores[rankP]
		fitness<-2-SelPress+(2*(SelPress-1)*(c(1:PopSize)-1)/(PopSize-1))
		
		#selection:stochastic uniform sampling 
		wheel1<-cumsum(fitness/sum(fitness))
		breaks<-runif(1)*1/PopSize
		breaks<-c(breaks,breaks+((1:(PopSize-1)))/PopSize)
		sel<-rep(1,PopSize)
		
		for(i in 1:length(breaks)){
			sel[i]<-which(wheel1>breaks[i])[1]
			}
			
		#intermediate generation
		Pop2<-Pop[sel,]
		PSize2<-dim(Pop2)[1]
		PSize3<-PopSize-elitism
		
		#Recombination: uniform: each bit has a .5 proba of being inherited from each parent
		mates<-cbind(ceiling(runif(PSize3)*PSize2),ceiling(runif(PSize3)*PSize2))
		
		#This holds the probability, for each bit, to be inherited from parent 1 (if TRUE) or 2 (if FALSE)
		InhBit<-matrix(runif((PSize3*bLength)),nrow=PSize3,ncol=bLength)
		InhBit<-InhBit < 0.5
		
		#Try one point crossover
		#xover<-ceiling(runif(PSize3)*(bLength-1))
		#indices<-matrix(1:bLength,nrow=PSize3,ncol=bLength,byrow=TRUE)
		#InhBit<-matrix(rep(FALSE,PSize3*bLength),nrow=PSize3,ncol=bLength)
		#for(i in 1:PSize3){
		#	InhBit[i,]<-indices[i,]<xover[i]
		#	}
		#
		
		Pop3par1<-Pop2[mates[,1],]
		Pop3par2<-Pop2[mates[,2],]
		Pop3<-Pop3par2
		Pop3[InhBit]<-Pop3par1[InhBit]
		
		#Mutation
		MutProba<-matrix(runif((PSize3*bLength)),nrow=PSize3,ncol=bLength)
		MutProba<-(MutProba < (Pmutation/bLength))
		Pop3[MutProba]<-1-Pop3[MutProba]
		
		#Compute stats
		t<-c(t,Sys.time())
		g<-g+1
		thisGenBest<-scores[length(scores)]
		thisGenBestBit<-Pop[length(scores),]
		
		if(is.na(thisGenBest)){
		
			thisGenBest<-min(scores, na.rm=TRUE)
			thisGenBestBit<-Pop[which(scores == thisGenBest)[1],]
			
			}
			
		if(thisGenBest < bestobj){
		
			bestobj<-thisGenBest
			bestbit<-thisGenBestBit
			stallGen<-0
			
			}else{
			
				stallGen<-stallGen+1
				}
				
		resThisGen<-c(
			g,
			bestobj,
			toString(bestbit),
			stallGen,
			(mean(scores,na.rm=TRUE)),
			thisGenBest,
			toString(thisGenBestBit),
			as.numeric((t[length(t)]-t[length(t)-1]), units="secs"))	
			
		names(resThisGen)<-c(
			"Generation",
			"Best_score",
			"Best_bitString",
			"Stall_Generation",
			"Avg_Score_Gen",
			"Best_score_Gen",
			"Best_bit_Gen",
			"Iter_time")
			
		if(verbose) print(resThisGen)
		
		res<-rbind(res,resThisGen)
		
		#Check stopping criteria
		Criteria<-c((stallGen > StallGenMax),(as.numeric((t[length(t)]-t[1]), units="secs") > MaxTime),(g > maxGens))
		if(any(Criteria)) stop<-TRUE
		
		#Check for bitstrings that are within the tolerance of the best bitstring
		tolScore<-scores[length(scores)]*RelTol
		TolBs<-which(scores < scores[length(scores)]+tolScore)
		
		if(length(TolBs) > 0){
			PopTol<-rbind(PopTol,Pop[TolBs,])
			PopTolScores<-c(PopTolScores,scores[TolBs])
			}
			
		if(elitism > 0){
			Pop<-rbind(Pop3,Pop[(PopSize-elitism+1):PopSize,])
			}else{
				Pop<-Pop3
				}
		}
#end of the while loop

	PopTol<-PopTol[-1,]
	PopTolScores<-PopTolScores[-1]
	TolBs<-which(PopTolScores < scores[length(scores)]+tolScore)
	PopTol<-PopTol[TolBs,]
	PopTolScores<-PopTolScores[TolBs]
	PopTolT<-cbind(PopTol,PopTolScores)
	PopTolT<-unique(PopTolT,MARGIN=1)
	
	if(!is.null(dim(PopTolT))){ 
		PopTol<-PopTolT[,1:(dim(PopTolT)[2]-1)]
		PopTolScores<-PopTolT[,dim(PopTolT)[2]]
		}else{
			PopTol<-PopTolT[1:(length(PopTolT)-1)]
			PopTolScores<-PopTolT[length(PopTolT)]
			}
			
	res<-res[3:dim(res)[1],]	
	rownames(res)<-NULL
	return(list(
		bString=bestbit,
		Results=res,
		StringsTol=PopTol,
		StringsTolScores=PopTolScores))	
	}


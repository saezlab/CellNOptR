#
#  This file is part of the CNO software
#
#  Copyright (c) 2011-2012 - EMBL - European Bioinformatics Institute
#
#  File author(s): CNO developers (cno-dev@ebi.ac.uk)
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  CNO website: http://www.cellnopt.org
#
##############################################################################
# $Id: gaBinaryT1.R 3855 2013-07-29 10:36:59Z cokelaer $
gaBinaryT1<-function(
    CNOlist,
    model,
    initBstring=NULL,
    sizeFac=0.0001,
    NAFac=1,
    popSize=50,
    pMutation=0.5,
    maxTime=60,
    maxGens=500,
    stallGenMax=100,
    selPress=1.2,
    elitism=5,
    relTol=0.1,
    verbose=TRUE,
    priorBitString=NULL,
    timeIndex=2){

    # by default initial bit string is made of ones.
    if (is.null(initBstring)==TRUE){
        initBstring<-rep(1,length(model$reacID))
    }

    if ((class(CNOlist)=="CNOlist")==FALSE){
        CNOlist = CellNOptR::CNOlist(CNOlist)
    }

    # should be after CNOlist conversion
    if (timeIndex<2){ stop("timeIndex must be >=2")}
    if (timeIndex>length(CNOlist@timepoints)){ 
        stop(paste("timeIndex must be <= ", length(CNOlist@timepoints), sep=" "))
    }

    # ---- section related to T1  ----
    bLength<-length(initBstring)

    simList = prep4sim(model)
    indexList = indexFinder(CNOlist, model)

    Pop<-rbind(
        initBstring,
        round(matrix(runif(bLength*(popSize-1)), 
        nrow=(popSize-1),ncol=bLength))
	)
    # ---- section related to T1  end ----

    Pop <- addPriorKnowledge(Pop, priorBitString)

    bestbit<-Pop[1,]
    bestobj<-Inf
    stop<-FALSE
    g<-0
    stallGen<-0
    res<-rbind(
        c(g,bestobj,toString(bestbit),stallGen,Inf,Inf,toString(bestbit),0),
        c(g,bestobj,toString(bestbit),stallGen,Inf,Inf,toString(bestbit),0))
    colnames(res)<-c("Generation","Best_score","Best_bitString","Stall_Generation",
        "Avg_Score_Gen","Best_score_Gen","Best_bit_Gen","Iter_time")
    PopTol<-rep(NA,bLength)
    PopTolScores<-NA

    library(hash)
    scores2Hash = hash()

    #Function that produces the score for a specific bitstring
    getObj<-function(x){

        key = toString(.int2dec(x))
        if (has.key(key, scores2Hash)==TRUE){
            return(scores2Hash[[key]])
        } else {
            Score = computeScoreT1(CNOlist, model, x, simList, indexList,
                sizeFac, NAFac, timeIndex)
            if (length(scores2Hash)<1000){
                scores2Hash[[key]] =  Score
            }
        }

        return(Score)
    }

    #Loop
    t0<-Sys.time()
    t<-t0



    if (popSize*stallGenMax > 2**bLength){
        print("Given your input parameter, an exhaustive search will be faster...")
        # better to perform an exhaustive search
        stop = TRUE # stop criteria of the GA that need not to be run
        res = exhaustive(CNOlist, model, relTol=relTol, 
                   sizeFac=sizeFac, NAFac=NAFac, verbose=verbose)

        return(res)
    }


    # if bitstring has only 1 bit to optimize, enter this simple loop:
    # we should have an exhaustive optimisation as well for simple cases.
    if (bLength==1){
        # build a population made of 2 vectors: c(0) and c(1)
        Pop = matrix(c(0,1), nrow=2)
        scores<-apply(Pop,1,getObj, scoresHash=NULL)  
        rankP<-order(scores,decreasing=TRUE)
        # extract the best solution.
        iBest = rankP[2]
        return(list(
            bString=Pop[iBest,],
            bScore=scores[iBest],
            results=res,
            stringsTol=PopTol,
            stringsTolScores=PopTolScores))
    }

    while(!stop){

        #compute the scores
        #l1 = length(scores2Hash)
        scores <- apply(Pop, 1, getObj)
        #l2 = length(scores2Hash)
        #if (verbose == TRUE){print(paste("new scores:", l2-l1))}


        #Fitness assignment: ranking, linear
        rankP<-order(scores,decreasing=TRUE)
        Pop<-Pop[rankP,]
        scores<-scores[rankP]
        fitness<-2-selPress+(2*(selPress-1)*(c(1:popSize)-1)/(popSize-1))

        #selection:stochastic uniform sampling
        wheel1<-cumsum(fitness/sum(fitness))
        breaks<-runif(1)*1/popSize
        breaks<-c(breaks,breaks+((1:(popSize-1)))/popSize)
        sel<-rep(1,popSize)

        for(i in 1:length(breaks)){
            sel[i]<-which(wheel1>breaks[i])[1]
            }

        #intermediate generation
        Pop2<-Pop[sel,]
        PSize2<-dim(Pop2)[1]
        PSize3<-popSize-elitism

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
        #    InhBit[i,]<-indices[i,]<xover[i]
        #    }
        #

        Pop3par1<-Pop2[mates[,1],]
        Pop3par2<-Pop2[mates[,2],]
        Pop3<-Pop3par2
        Pop3[InhBit]<-Pop3par1[InhBit]

        #Mutation
        MutProba<-matrix(runif((PSize3*bLength)),nrow=PSize3,ncol=bLength)
        MutProba<-(MutProba < (pMutation/bLength))
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

        names(resThisGen)<-c("Generation","Best_score","Best_bitString","Stall_Generation",
            "Avg_Score_Gen","Best_score_Gen","Best_bit_Gen","Iter_time")

        if(verbose) {
            this = resThisGen
            this[[3]] = substring(this[[3]], 1, 80)
            this[[7]] = substring(this[[7]], 1, 80)
            print(this)
        }


        res<-rbind(res,resThisGen)

        #Check stopping criteria
        Criteria<-c((stallGen > stallGenMax),
            (as.numeric((t[length(t)]-t[1]), units="secs") > maxTime),
            (g > maxGens))
        if(any(Criteria)) stop<-TRUE

        #Check for bitstrings that are within the tolerance of the best bitstring
        tolScore<-scores[length(scores)]*relTol
        TolBs<-which(scores < scores[length(scores)]+tolScore)

        if(length(TolBs) > 0){
            PopTol<-rbind(PopTol,Pop[TolBs,])
            PopTolScores<-c(PopTolScores,scores[TolBs])
            }

        if(elitism > 0){
            Pop<-rbind(Pop3,Pop[(popSize-elitism+1):popSize,])
            }else{
                Pop<-Pop3
                }
        Pop <- addPriorKnowledge(Pop, priorBitString)
    }
    #end of the while loop

    PopTol<-as.matrix(PopTol[-1,])
    PopTolScores<-PopTolScores[-1]
    TolBs<-which(PopTolScores < scores[length(scores)]+tolScore)
    PopTol<-as.matrix(PopTol[TolBs,])
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
        bScore=bestobj,
        results=res,
        stringsTol=PopTol,
        stringsTolScores=PopTolScores))
    }


addPriorKnowledge <- function(pop, priorBitString){
    if (is.null(priorBitString) == TRUE){
        return(pop)
    }
    else{
        for (i in 1:dim(pop)[1]){
            pop[i,!is.na(priorBitString)] = priorBitString[!is.na(priorBitString)]
        }
    }
   return(pop)
}


.int2dec <- function(x){
    return(sum(x*2^(rev(seq(x))-1)))
}

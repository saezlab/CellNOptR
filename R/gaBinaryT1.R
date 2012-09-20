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
# $Id: gaBinaryT1.R 2455 2012-09-20 09:20:00Z cokelaer $
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
    maxSizeHashTable=5000){

    # by default initial bit string is made of ones.
    if (is.null(initBstring)==TRUE){
        initBstring<-rep(1,length(model$reacID))
    }

    if ((class(CNOlist)=="CNOlist")==FALSE){
        CNOlist = CellNOptR::CNOlist(CNOlist)
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
    obj<-rep(0,popSize)
    g<-0
    stallGen<-0
    res<-rbind(
        c(g,bestobj,toString(bestbit),stallGen,Inf,Inf,toString(bestbit),0),
        c(g,bestobj,toString(bestbit),stallGen,Inf,Inf,toString(bestbit),0))
    colnames(res)<-c("Generation","Best_score","Best_bitString","Stall_Generation",
        "Avg_Score_Gen","Best_score_Gen","Best_bit_Gen","Iter_time")
    PopTol<-rep(NA,bLength)
    PopTolScores<-NA

    #Function that produces the score for a specific bitstring
    getObj<-function(x, scoresHash=NULL){

        bitString<-x

        # the hash table is used to speed up code. gain is guaranteed to be at least equal to elitism/popsize
        if (is.null(scoresHash)==FALSE){
            thisScore <- scoresHash[rownames(scoresHash) == paste(unlist(x), collapse=","),1]
             if (length(thisScore) != 0){
                 return(thisScore)
            } # otherwise let us keep going
        }

        Score = computeScoreT1(CNOlist, model, bitString, simList, indexList,
			sizeFac, NAFac)

        return(Score)
    }

    #Loop
    t0<-Sys.time()
    t<-t0

    # used by the scores hashTable.
    scoresHash <- data.frame()
    # if you do want the hastable, uncomment the following line.
    #scoresHash = NULL

    while(!stop){

        #compute the scores
        scores<-apply(Pop,1,getObj, scoresHash=scoresHash)

        # fill the hash table to speed up code
        scoresHash<-fillHashTable(scoresHash, scores, Pop, maxSizeHashTable)

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

        if(verbose) print(resThisGen)

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

# simple function to shift a data.frame
shift <- function(d, k) rbind( tail(d,k), head(d,-k), deparse.level = 0 )



fillHashTable <-function(scoresHash, scores, Pop, maxSizeHashTable=5000)
{
    # if not a data.frame, just return NULL
    if (is.null(scoresHash)==TRUE){ return(NULL)}

    popSize = dim(Pop)[1]
    for (i in 1:dim(Pop)[1]){
        thisScore <- scoresHash[rownames(scoresHash) == paste(unlist(Pop[i,]), collapse=","), 1]
        # if not already stored, store the score and corresponding bitstring
        if (length(thisScore) == 0){
            # compute a new score
            thisScore <- scores[i]
            # rename the row (latest one) of the newly added score
            if (dim(scoresHash)[1]<maxSizeHashTable){
                scoresHash <- rbind(scoresHash, c(thisScore, 1))
                j = dim(scoresHash)[1]
                row.names(scoresHash)[j] = paste(unlist(Pop[i,]), collapse=",")
            }
            else{
                # shift by -1 so that first element put in the queue (that
                # we get rid of) is last
                #indices = c(1:maxSizeHashTable-popSize*2)
                #scoresHash[indices, ] = scoresHash[order(scoresHash[indices,2], decreasing=TRUE), ]
                scoresHash = shift(scoresHash, -1)


                #scoresHash = shift(scoresHash, -1)
                # overwrite last element with the latest score and bitstring
                scoresHash[maxSizeHashTable,] = c(thisScore, 1)
                row.names(scoresHash)[maxSizeHashTable] =
                    paste(unlist(Pop[i,]), collapse=",")
            }
         }
         else {
             count = scoresHash[rownames(scoresHash) == paste(unlist(Pop[i,]), collapse=","), 2]
             scoresHash[rownames(scoresHash) == paste(unlist(Pop[i,]), collapse=","), 2] = count+1
         }
    }
    return(scoresHash)
}

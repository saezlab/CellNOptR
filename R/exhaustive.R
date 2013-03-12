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
# $Id$
exhaustive<-function(
    CNOlist,
    model,
    shuffle=FALSE,
    Nmax=NULL,
    verbose=TRUE,
    sizeFac = 0.0001,
    NAFac = 1,
    relTol=0.1,
    timeIndex=2){

    if ((class(CNOlist)=="CNOlist")==FALSE){
        CNOlist = CellNOptR::CNOlist(CNOlist)
    }

    # should be after CNOlist conversion
    if (timeIndex<2){ stop("timeIndex must be >=2")}
    if (timeIndex>length(CNOlist@timepoints)){ 
        stop(paste("timeIndex must be <= ", length(CNOlist@timepoints), sep=" "))
    }

    bLength = length(model$reacID)
    if (bLength>20){
        print(paste("You will compute ", 2**bLength, "iterations. It may take a
            while..."), sep="")
    }
    if (bLength == 0){
        stop("Nothing to optimise...bitstring length is zero.")
    }
    # boolean case but could be extended easily to fuzzy by changing
    # rep(list(0,1)) to rep(list(0, nTF))
    bitstrings = expand.grid(rep(list(0:1), bLength)) 

    # do we want to shuffle the order ?
    if (shuffle==TRUE){
        bitstrings = bitstrings[sample.int(nrow(bitstrings)),]
    }

    # compute these variables once for all
    simList = prep4sim(model)
    indexList = indexFinder(CNOlist, model)

    bestScore = 1e6
    bitstring = bitstrings[1]
    N = 2**bLength

    # if the user request a Nmax value, let us use it
    if (is.null(Nmax)==FALSE){
        if (Nmax<=0){stop("Nmax must be positive stricly")}
        N = min(Nmax, N)
    } else{
        Nmax=N
    }

    # NA are removed later on
    PopTol<-rep(NA,bLength)
    PopTolScores<-NA 

    all_scores = c()
    t1 = Sys.time()
    for (x in seq(1:Nmax)){

        bitstring = as.double(bitstrings[x,])
        #if (sum(bitstring)==0){next}
        score = computeScoreT1(CNOlist, model, bitstring, simList, indexList,
			sizeFac, NAFac, timeIndex)
        all_scores[x] = score
        changed = FALSE
        if (score <= bestScore){
            if (verbose==TRUE){
                print(paste("--Found a new best score=", score, " at iteration ", x, sep=""))
                print(bitstring)
            }
            bestScore = score
            bestbit = bitstring
            changed = TRUE
        } 

        if (verbose == TRUE){
            if ((x%%1000)==0){
                print(paste(x,"/",N,sep=""))
            }
        }

        # Check if the current score/bitstring is within tolerance of the best
        # score
        if( score <=(1+relTol)*bestScore){

            PopTol<-rbind(PopTol, bitstring)
            PopTolScores<-c(PopTolScores, score)
        }

        # if the best score changed, there are maybe scores to be removed now
        if (changed==TRUE){
            indices = PopTolScores <= bestScore* (1+relTol)
            PopTol = PopTol[indices,]
            PopTolScores = PopTolScores[indices]
        }

    }
    t2 = Sys.time()
    if (verbose==TRUE){
        print(t2-t1)    
    }

    # remove names of the matrix row and column
    PopTol = matrix(PopTol, ncol=bLength, 
        dimnames=list(rep("", dim(PopTol)[1]),rep("",bLength)))

    # remove first NA element. The NA element has the good property that if
    # there is only one score found within the tol, there are 2 elements hence
    # we are still dealing with a matrix. So, we remove the NA only at the end.

    # removing the two rows.
    PopTol = matrix(PopTol[-1,], ncol=bLength)
    PopTolScores = matrix(PopTolScores[-1])

    # todo
    res = NULL

    return(list(
        bString=bestbit,
        bScore=bestScore,
        results=res,
        stringsTol=PopTol,
        stringsTolScores=PopTolScores,
        scores=all_scores))
    }



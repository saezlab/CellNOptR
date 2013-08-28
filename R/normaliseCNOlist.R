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
# $Id: normaliseCNOlist.R 3676 2013-06-05 12:27:59Z cokelaer $
normaliseCNOlist <- function(
    CNOlist,
    EC50Data=0.5,
    HillCoef=2,
    EC50Noise=0.,
    detection=0,
    saturation=Inf,
    changeTh=0,
    norm2TorCtrl=NULL,
    mode="time", options=list(rescale_negative=T),
    verbose=FALSE){


    if ((class(CNOlist)=="CNOlist")==FALSE){
         CNOlist = CellNOptR::CNOlist(CNOlist)
     }

    if (is.null(norm2TorCtrl)==FALSE){
        stop("'norm2TorCtrl' argument is deprecated. Please use 'mode=' instead")
    }

    #Check that the parameters have the right format
    if(class(EC50Data) != "numeric" | length(EC50Data) != 1)
        warning("The parameter 'EC50Data' should be a single numeric value")
    if(class(EC50Noise) != "numeric" | length(EC50Noise) != 1)
        warning("The parameter 'EC50Noise' should be a single numeric value")
    if((class(HillCoef) != "numeric" & class(HillCoef) != "integer") | length(HillCoef) != 1)
        stop("The parameter 'HillCoef' should be a single numeric value")
    if((class(detection) != "numeric" & class(detection) != "integer")  | length(detection) != 1)
        stop("The parameter 'detection' should be a single numeric value")
    if((class(changeTh) != "numeric" & class(changeTh) != "integer") | length(changeTh) != 1)
        stop("The parameter 'changeTh' should be a single numeric value")
    if((class(saturation) != "numeric" & class(saturation) != "integer") | length(saturation) != 1)
        stop("The parameter 'Saturation' should be a single numeric value")
    if(mode %in% c("time","ctrl", "raw") == FALSE)
        stop("The parameter 'mode' should be either 'time' or 'ctrl' or 'raw'")
    

    if (verbose==TRUE){
        cat("Normalisation mode is: ", mode, "\n", sep="")
    }
    # Check which values are out of the dynamic range of the machine: create a list
    # of matrices, one matrix for each time point (i.e. each element of the field
    # CNOlist$valueSignals) filled with FALSE if the measurement is in the dynamic
    # range, and true if it isn't (i.e. it should be replaced by NA at the end) i.e.
    # if Detection=0 and saturation=Inf, then all the matrices in NaNsList should be
    # filled with 0s this function also tags as NAs all the field that don't have a
    # value (should be imported as NA)

    NaNsList <- CNOlist@signals

    DynamicRange <- function(x){
        if (verbose==TRUE){
            data_sat = x[which(x > saturation )]
            cat("found ", length(data_sat), " data points above saturation (", saturation, ").\n", sep="")
        }
        x[which(x > saturation )] <- Inf

        if (verbose==TRUE){
            data = x[which(x < detection)]
            cat("found ", length(data), " data points below detection (",detection, ").\n", sep="")
        }
        x[which(x < detection )] <- Inf

        x[which(is.na(x))] <- Inf
        x[is.finite(x)] <- FALSE
        x[is.infinite(x)] <- TRUE
        return(x)
        }

    NaNsList <- lapply(NaNsList, DynamicRange )
    # Calculate relative change

    # the following list will contain true for the fold changes that are negatives
    # and significant,and false for the other ones.
    negList <- CNOlist@signals
    negList <- lapply(negList,
                    function(x) x <- matrix(data=FALSE,
                    nrow=dim(x)[1],ncol=dim(x)[2]) 
                )

    # the following list will contain true for the fold changes that are positive
    # and significant, and false for the other ones.
    posList <- CNOlist@signals
    posList <- lapply(posList,
                    function(x) x <- matrix(data=FALSE,
                    nrow=dim(x)[1],ncol=dim(x)[2]) 
                )

    # the following matrix will contain the actual fold change values
    FoldChangeList <- CNOlist@signals

    # if mode="raw" data may be negative
    if(mode == "raw"){
        for(i in 2:length(FoldChangeList)){
            negList[[i]][which((CNOlist@signals[[i]] - CNOlist@signals[[1]]) < (-1*changeTh))] <- TRUE
            posList[[i]][which((CNOlist@signals[[i]] - CNOlist@signals[[1]]) > changeTh)] <- TRUE
        }
    }

    # if mode="time" then the relative change is simply matrix t1 - matrix
    # t0 (ie each measurement is compared to the exact same condition at time 0).

    if(mode == "time"){
        for(i in 2:length(FoldChangeList)){
            FoldChangeList[[i]] <- abs(CNOlist@signals[[i]] - CNOlist@signals[[1]])/CNOlist@signals[[1]]


            negList[[i]][which((CNOlist@signals[[i]] - CNOlist@signals[[1]]) < (-1*changeTh))] <- TRUE
            posList[[i]][which((CNOlist@signals[[i]] - CNOlist@signals[[1]]) > changeTh)] <- TRUE
        }

        FoldChangeList[[1]] <- matrix(
            0,
            ncol=dim(FoldChangeList[[1]])[2],
            nrow=dim(FoldChangeList[[1]])[1])

        }

    # else
    if(mode == "ctrl"){

    # if mode="ctrl" then the relative change is computed relative to the
    # ctrl at the same time the ctrl is the case without stimuli, but with the same
    # inhibitors.   In our case this still means that at t0 we are going to have
    # zero everywhere since only two measurements were made: with and without the
    # inhibitor(s) and these measurements have been copied across corresponding
    # position this last bit makes sense because we assume that the inhibitors are
    # already present at time 0 when we add the stimuli to find the right row to
    # normalise, we look for a row in the valueStimuli that has only 0's but where
    # the corresponding row in the inhibitor matrix has the same status as the row
    # that we are trying to normalise.
        for(i in 2:length(FoldChangeList)){
            for(n in 1:dim(FoldChangeList[[i]])[1]){
                #First I treat the rows that are not ctrls
                if(sum(CNOlist@stimuli[n,]) != min(apply(CNOlist@stimuli,MARGIN=1,sum))){
                    ctrlRow <- intersect(
                        which(
                            apply(CNOlist@stimuli,MARGIN=1,sum) ==
                                min(apply(CNOlist@stimuli,MARGIN=1,sum))),
                        which(
                            apply(CNOlist@inhibitors,MARGIN=1,
                                function(x) all(x == CNOlist@inhibitors[n,]))))
                    FoldChangeList[[i]][n,] <- abs(CNOlist@signals[[i]][n,] - CNOlist@signals[[i]][ctrlRow,])/CNOlist@signals[[i]][ctrlRow,]
                    negList[[i]][n,which((CNOlist@signals[[i]][n,] - CNOlist@signals[[i]][ctrlRow,]) < (-1*changeTh) )] <- TRUE
                    posList[[i]][n,which((CNOlist@signals[[i]][n,] - CNOlist@signals[[i]][ctrlRow,]) > changeTh )] <- TRUE

                    }else{

                        #Then I set to 0 all the rows that are ctrls
                        FoldChangeList[[i]][n,] <- rep(0,dim(FoldChangeList[[i]])[2])
                    }

                }
            }

            FoldChangeList[[1]] <- matrix(
                0,
                ncol=dim(FoldChangeList[[1]])[2],
                nrow=dim(FoldChangeList[[1]])[1])

        }
    # Now I compute the penalty for being noisy, which is calculated for each
    # measurement as the measurement divided by the max measurement across all
    # conditions and time for that particular signal (excluding values out of the
    # dynamic range).

    #1.This small function computes the max across all measurements for each signal,
    #excluding the values out of the dynamic range
    SignalMax <- function(signals,NaNsList){

        for(i in 1:length(signals)){
            signals[[i]][which(NaNsList[[i]] == 1)] <- 0
            }

        for(i in 2:length(signals)){
            signals[[i]] <- rbind(signals[[i-1]],signals[[i]])
            }

        signals <- signals[[length(signals)]]
        if (mode=="raw"){ # in raw mode, negative values are possible, need to take abs
            signalABS <- apply(signals,MARGIN=2,abs)
            signalMax <- apply(signalABS,MARGIN=2,max)
        }else{
            signalMax <- apply(signals,MARGIN=2,max)
        }
        return(signalMax)
        }

    signalMax <- SignalMax(signals=CNOlist@signals,NaNsList=NaNsList)
    # 2.This bit takes a list of matrices and goes into each matrix
    # and divides each row by the vector containing the max values
    PenaltyList <- CNOlist@signals
    PenaltyList <- lapply(
        PenaltyList,
        function(x){
            x <- apply(x,MARGIN=1,function(y){y <- y/signalMax});
            return(t(x))
        })

    #3.Now I can make this list of values go through the saturation function
    if (mode=="raw"){
        # in raw mode, negative values are possible. Just take the absolute values
        SatPenalty <- lapply(PenaltyList,abs)
    }else{
        SatPenalty <- PenaltyList
     }
    SatPenalty <- lapply(SatPenalty,function(x) {x/(EC50Noise+x)} )

    # Now I make the data go through the Hill function
    HillData <- FoldChangeList
    HillData <- lapply(
        HillData,
        function(x) {x^HillCoef/((EC50Data^HillCoef)+(x^HillCoef))} )

    # Now I can multiply HillData and SatPenalty, matrix by matrix and element by
    # element.
    NormData <- HillData
    for(i in 1:length(NormData)){
        NormData[[i]] <- HillData[[i]]*SatPenalty[[i]]
        NormData[[i]][which(negList[[i]] == 1)] <- NormData[[i]][which(negList[[i]] == 1)]*-1
        if(mode != "raw"){ # in raw mode, neg+pos can leaad to zero. time 0 is also zero.
            NormData[[i]][which((negList[[i]] + posList[[i]]) == 0)] <- 0
        }
    }

    # Now let us set to NaN all the values that have been tagged in the beginning as
    # out of the dynamic range.
    for(i in 1:length(NormData)){
        NormData[[i]][which(NaNsList[[i]] == 1)] <- NaN
    }

    # rescale the columns that have negative values. T0 is untouched.
    if (options$rescale_negative == T && 1==0){ # this loop does not work. added 1==0 temporary so that we do not enter in the loop
        for (x in colnames(NormData[[1]])){
            for (i in 2:length(NormData)){ #there is always a Time 0 so were are guarantee to have at least 2 time points
                m = min(NormData[[i]][,x], na.rm=T)
                M = max(NormData[[i]][,x], na.rm=T)
                if (m<0 && m!=M ){
                    NormData[[i]][,x] = (NormData[[i]][,x]-m)/ (M-m)
                }
            }
        }
    }

    # The following code is ONE way of dealing with negative values
    # !!! each experiment/specy is treated separately.
    # 0,-0.5,-1 will be rescaled in 1,0.5,0
    # !!! 0,-0.01,0.5 is also rescaled because there is 1 negative value
    if (options$rescale_negative == T){
        # prepare a more convenient data set to manipulate

        # rescale each column independently
        for (x in colnames(NormData[[1]])){
            # extract only data related to the specy x
            c = NormData[[1]][,x]
            for (i in 2:length(NormData)){
                c = rbind(c, NormData[[i]][,x])
            }
            # get the min and max over time (apply) for each experiment
            min_vector = apply(c, 2, min, na.rm=T)
            max_vector = apply(c, 2, max, na.rm=T)

            # rescale negative values if needed for each experiment
            for (i in seq_along(min_vector)){
                m = min_vector[i]
                M = max_vector[i]
                if (m<0 && m!=M){
                    c[,i] = (c[,i]-m)/(M-m)
                    # Finally, save the new values in NormData
                    for (itime in 1:length(NormData)){
                        NormData[[itime]][i,x] = c[itime,i]
                    }
                }
            }

        }
    } 


    CNOlist@signals <- NormData

    # Variance should be zero
    for (i in seq_along(CNOlist@timepoints)){
        CNOlist@variances[[i]] = CNOlist@variances[[i]] * 0
    }

    return(CNOlist)
}


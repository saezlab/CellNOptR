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
# $Id: normaliseCNOlist.R 1586 2012-06-26 14:59:24Z cokelaer $
normaliseCNOlist<-function(
    CNOlist,
    EC50Data=0.5,
    HillCoef=2,
    EC50Noise=0.1,
    detection=0,
    saturation=Inf,
    changeTh=0,
    norm2TorCtrl="time"){

#check that CNOlist is a CNOlist
    if(!is.list(CNOlist)){
        stop("This function expects as input a CNOlist as output by makeCNOlist")
        }
    if(all(names(CNOlist) != c(
        "namesCues",
        "namesStimuli",
        "namesInhibitors",
        "namesSignals",
        "timeSignals",
        "valueCues",
        "valueInhibitors",
        "valueStimuli",
        "valueSignals"))){
        stop("This function expects as input a CNOlist as output by makeCNOlist")
        }

#Check that the parameters have the right format
    if(class(EC50Data) != "numeric" | length(EC50Data) != 1)
        warning("The parameter 'EC50Data' should be a single numeric value")
    if(class(EC50Noise) != "numeric" | length(EC50Noise) != 1)
        warning("The parameter 'EC50Noise' should be a single numeric value")
    if((class(HillCoef) != "numeric" & class(HillCoef) != "integer") | length(HillCoef) != 1)
        warning("The parameter 'HillCoef' should be a single numeric value")
    if((class(detection) != "numeric" & class(detection) != "integer") | length(detection) != 1)
        warning("The parameter 'detection' should be a single numeric value")
    if((class(changeTh) != "numeric" & class(changeTh) != "integer") | length(changeTh) != 1)
        warning("The parameter 'changeTh' should be a single numeric value")
    if((class(saturation) != "numeric" & class(saturation) != "integer") | length(saturation) != 1)
        warning("The parameter 'Saturation' should be a single numeric value")
    if(length(norm2TorCtrl) != 1 | length(match(norm2TorCtrl, c("time","ctrl"))) != 1)
        warning("The parameter 'norm2TorCtrl' should be either 'time' or 'ctrl'")


#Check which values are out of the dynamic range of the machine: create a list of matrices,
#one matrix for each time point (i.e. each element of the field CNOlist$valueSignals)
#filled with FALSE if the measurement is in the dynamic range, and true if it isn't
#(i.e. it should be replaced by NA at the end)
#i.e. if Detection=0 and saturation=Inf, then all the matrices in NaNsList should be filled with 0s
#this function also tags as NAs all the field that don't have a value (should be imported as NA)

    NaNsList<-CNOlist$valueSignals

    DynamicRange<-function(x){
        x[which(x > saturation )]<-Inf
        x[which(x < detection )]<-Inf
        x[which(is.na(x))]<-Inf
        x[is.finite(x)]<-FALSE
        x[is.infinite(x)]<-TRUE
        return(x)
        }

    NaNsList<-lapply(NaNsList, DynamicRange )

#Calculate relative change

#the following list will contain true for the fold changes that are negatives and significant,
#and false for the other ones
    negList<-CNOlist$valueSignals
    negList<-lapply(negList,function(x) x<-matrix(data=FALSE,nrow=dim(x)[1],ncol=dim(x)[2]) )

#the following list will contain true for the fold changes that are positive and significant,
#and false for the other ones
    posList<-CNOlist$valueSignals
    posList<-lapply(posList,function(x) x<-matrix(data=FALSE,nrow=dim(x)[1],ncol=dim(x)[2]) )

#the following matrix will contain the actual fold change values
    FoldChangeList<-CNOlist$valueSignals

#if norm2TorCtrl="time" then the relative change is simply matrix t1 - matrix t0
#(ie each measurement is compared to the exact same condition at time 0)
    if(norm2TorCtrl == "time"){

        for(i in 2:length(FoldChangeList)){
            FoldChangeList[[i]]<-abs(CNOlist$valueSignals[[i]] - CNOlist$valueSignals[[1]])/CNOlist$valueSignals[[1]]
            negList[[i]][which((CNOlist$valueSignals[[i]] - CNOlist$valueSignals[[1]]) < (-1*changeTh))]<-TRUE
            posList[[i]][which((CNOlist$valueSignals[[i]] - CNOlist$valueSignals[[1]]) > changeTh)]<-TRUE
            }

        FoldChangeList[[1]]<-matrix(
            0,
            ncol=dim(FoldChangeList[[1]])[2],
            nrow=dim(FoldChangeList[[1]])[1])

        }else{

#if norm2TorCtrl="ctrl" then the relative change is computed relative to the ctrl at the same time
#the ctrl is the case without stimuli, but with the same inhibitors.
#In our case this still means that at t0 we are going to have zero everywhere since
#only two measurements were made: with and without the inhibitor(s)
#and these measurements have been copied across corresponding position
#this last bit makes sense because we assume that the inhibitors are already present at time 0
#when we add the stimuli
#to find the right row to normalise, we look for a row in the valueStimuli that has only 0's but
#where the corresponding row in the inhibitor matrix has the same status as the row that we are trying to normalise

        for(i in 2:length(FoldChangeList)){

            for(n in 1:dim(FoldChangeList[[i]])[1]){
#First I treat the rows that are not ctrls

                    if(sum(CNOlist$valueStimuli[n,]) != min(apply(CNOlist$valueStimuli,MARGIN=1,sum))){
                        ctrlRow<-intersect(
                            which(
                                apply(CNOlist$valueStimuli,MARGIN=1,sum) ==
                                    min(apply(CNOlist$valueStimuli,MARGIN=1,sum))),
                            which(
                                apply(CNOlist$valueInhibitors,MARGIN=1,
                                    function(x) all(x == CNOlist$valueInhibitors[n,]))))
                        FoldChangeList[[i]][n,]<-abs(CNOlist$valueSignals[[i]][n,] - CNOlist$valueSignals[[i]][ctrlRow,])/CNOlist$valueSignals[[i]][ctrlRow,]
                        negList[[i]][n,which((CNOlist$valueSignals[[i]][n,] - CNOlist$valueSignals[[i]][ctrlRow,]) < (-1*changeTh) )]<-TRUE
                        posList[[i]][n,which((CNOlist$valueSignals[[i]][n,] - CNOlist$valueSignals[[i]][ctrlRow,]) > changeTh )]<-TRUE

                        }else{

#Then I set to 0 all the rows that are ctrls
                            FoldChangeList[[i]][n,]<-rep(0,dim(FoldChangeList[[i]])[2])
                            }

                    }
                }

            FoldChangeList[[1]]<-matrix(
                0,
                ncol=dim(FoldChangeList[[1]])[2],
                nrow=dim(FoldChangeList[[1]])[1])

            }
#Now I compute the penalty for being noisy, which is calculated for each measurement
#as the measurement divided by the max measurement across all conditions and time
#for that particular signal (excluding values out of the dynamic range)

    #1.This small function computes the max across all measurements for each signal,
    #excluding the values out of the dynamic range
    SignalMax<-function(signals,NaNsList){

        for(i in 1:length(signals)){
            signals[[i]][which(NaNsList[[i]] == 1)]<-0
            }

        for(i in 2:length(signals)){
            signals[[i]]<-rbind(signals[[i-1]],signals[[i]])
            }

        signals<-signals[[length(signals)]]
        signalMax<-apply(signals,MARGIN=2,max)
        return(signalMax)
        }

    signalMax<-SignalMax(signals=CNOlist$valueSignals,NaNsList=NaNsList)

    #2.This bit takes a list of matrices and goes into each matrix
    #and divides each row by the vector containing the max values
    PenaltyList<-CNOlist$valueSignals
    PenaltyList<-lapply(
        PenaltyList,
        function(x){x<-apply(x,MARGIN=1,function(y){y<-y/signalMax});return(t(x))})

    #3.Now I can make this list of values go through the saturation function
    SatPenalty<-PenaltyList
    SatPenalty<-lapply(SatPenalty,function(x) {x/(EC50Noise+x)} )

#Now I make the data go through the Hill function
    HillData<-FoldChangeList
    HillData<-lapply(
        HillData,
        function(x) {x^HillCoef/((EC50Data^HillCoef)+(x^HillCoef))} )

#Now I can multiply HillData and SatPenalty, matrix by matrix and element by element
    NormData<-HillData

    for(i in 1:length(NormData)){
        NormData[[i]]<-HillData[[i]]*SatPenalty[[i]]
        NormData[[i]][which(negList[[i]] == 1)]<-NormData[[i]][which(negList[[i]] == 1)]*-1
        NormData[[i]][which((negList[[i]] + posList[[i]]) == 0)]<-0
        }

#Now I will set to NaN all the values that have been tagged in the beginning
#as out of the dynamic range

    for(i in 1:length(NormData)){
        NormData[[i]][which(NaNsList[[i]] == 1)]<-NaN
        }

    CNOlist$valueSignals<-NormData
    return(CNOlist)

    }


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
# $Id: getFit.R 2267 2012-08-30 15:31:54Z cokelaer $
getFit<-function(
    simResults,
    CNOlist,
    model,
    indexList,
    timePoint=c("t1","t2"),
    sizeFac=0.0001,
    NAFac=1,
    nInTot,
    simResultsT0=NA
    ){

    if ((class(CNOlist)=="CNOlist")==FALSE){
         CNOlist = CellNOptR::CNOlist(CNOlist)
     }


    simResults<-simResults[,indexList$signals]


    # for back compatibility, timePoint ca be "t1" or "t2" but developers should
    # use an integer.
    if(timePoint == "t1"){
        tPt<-2
    }
    else{
        if(timePoint == "t2"){
            tPt<-3
        }
        else{
            tPt<-timePoint
        }
    }


    # if t0 is provided and we are interested in t1
    # then  score is based on t1 but also t0
    if (tPt == 2 && is.na(simResultsT0)==FALSE){
        Diff0<-simResultsT0[,indexList$signals]-CNOlist@signals[[1]]
        Diff<-simResults-CNOlist@signals[[tPt]]
        r0<-Diff0^2
        r<-Diff^2
        r <- rbind(r0, r) # we can concatenate because it's matricial computation.
        deviationPen<-sum(r[!is.na(r)])/2
    }# otherwise, no need to take to into account
    else{
        Diff<-simResults-CNOlist@signals[[tPt]]
        r<-Diff^2
        deviationPen<-sum(r[!is.na(r)])
    }


    NAPen<-NAFac*length(which(is.na(simResults)))

    nDataPts<-dim(CNOlist@signals[[tPt]])[1]*dim(CNOlist@signals[[tPt]])[2]

    nInputs<-length(which(model$interMat == -1))

    # nInTot: number of inputs of expanded model
    # nInputs: number of inputs of cut model
    sizePen<-(nDataPts*sizeFac*nInputs)/nInTot

    score<-deviationPen+NAPen+sizePen

    return(score)

    }


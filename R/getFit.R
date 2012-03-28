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
# $Id: getFit.R 804 2012-03-22 16:56:26Z cokelaer $
getFit<-function(
	SimResults,
	CNOlist,
	Model,
	indexList,
	timePoint=c("t1","t2"),
	sizeFac=0.0001,
	NAFac=1,
	nInTot, 
    SimResultsT0=NA
    ){
	
	SimResults<-SimResults[,indexList$signals]
	
	if(timePoint == "t1") tPt<-2
	if(timePoint == "t2") tPt<-3

    # if t0 is provided and we are interested in t1
    # then  score is based on t1 but also t0
    if (tPt == 2 && is.na(SimResultsT0)==FALSE){
        Diff0<-SimResultsT0[,indexList$signals]-CNOlist$valueSignals$t0
        Diff<-SimResults-CNOlist$valueSignals[[tPt]]
    	r0<-Diff0^2
	    r<-Diff^2
        r <- rbind(r0, r) # we can concatenate because it's matricial computation.
	    deviationPen<-sum(r[!is.na(r)])/2
    }# otherwise, no need to take to into account
    else{
        Diff<-SimResults-CNOlist$valueSignals[[tPt]]
    	r<-Diff^2
    	deviationPen<-sum(r[!is.na(r)])
    }
    
	
	NAPen<-NAFac*length(which(is.na(SimResults)))
	
	nDataPts<-dim(CNOlist$valueSignals[[tPt]])[1]*dim(CNOlist$valueSignals[[tPt]])[2]
	
	nInputs<-length(which(Model$interMat == -1))
	
	# nInTot: number of inputs of expanded model
	# nInputs: number of inputs of cut model
	sizePen<-(nDataPts*sizeFac*nInputs)/nInTot
	
	score<-deviationPen+NAPen+sizePen
	
	return(score)
	
	}


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
# $Id: residualError.R 595 2012-02-22 17:21:47Z cokelaer $
residualError<-function(CNOlist){

#check that CNOlist is a CNOlist
	if(!is.list(CNOlist)){
		stop("This function expects as input a CNOlist as output by makeCNOlist or normaliseCNOlist")
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
		
	resErr<-rep(NA,3)
	names(resErr)<-c("t1","t2","t1andt2")
	Diff1<-round(CNOlist$valueSignals[[2]])-CNOlist$valueSignals[[2]]
	resErr[1]<-sum(Diff1^2,na.rm=TRUE)
	
	if(length(CNOlist$valueSignals) == 3){
		Diff2<-round(CNOlist$valueSignals[[3]])-CNOlist$valueSignals[[3]]
		resErr[2]<-sum(Diff2^2,na.rm=TRUE)
		resErr[3]<-resErr[1]+resErr[2]
		}
		
	if(length(CNOlist$valueSignals) > 3){
		warning("This version of the software only handles 2 time points")
		}
		
	return(resErr)	
	}


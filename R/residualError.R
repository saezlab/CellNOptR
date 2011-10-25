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


getFit<-function(
	SimResults,
	CNOlist,
	Model,
	indexList,
	timePoint=c("t1","t2"),
	sizeFac=0.0001,
	NAFac=1){
	
	SimResults<-SimResults[,indexList$signals]
	
	if(timePoint == "t1") tPt<-2
	if(timePoint == "t2") tPt<-3
	
	Diff<-SimResults-CNOlist$valueSignals[[tPt]]
	r<-Diff^2
	
	deviationPen<-sum(r[!is.na(r)])
	
	NAPen<-NAFac*length(which(is.na(SimResults)))
	
	nDataPts<-dim(CNOlist$valueSignals[[tPt]])[1]*dim(CNOlist$valueSignals[[tPt]])[2]
	
	nReac<-length(Model$reacID)
	
	nInputs<-length(which(Model$interMat == -1))
	
	sizePen<-(nDataPts*sizeFac*nInputs)/nReac
	
	score<-deviationPen+NAPen+sizePen
	
	return(score)
	
	}


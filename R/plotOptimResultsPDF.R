plotOptimResultsPDF<-function(
	SimResults=SimResults,
	expResults=expResults,
	times=times,
	namesCues=namesCues,
	namesSignals=namesSignals,
	valueCues=valueCues,
	fileName){
	
	if(sum(dim(SimResults[[1]])) < 20){
	
		pdf(file=fileName,width=14,height=7)
		
		}else{
		
			pdf(file=fileName,width=21,height=10)
			
			}
			
	plotOptimResults(
		SimResults=SimResults,
		expResults=expResults,
		times=times,
		namesCues=namesCues,
		namesSignals=namesSignals, 
		valueCues=valueCues)
	
	dev.off()
	}


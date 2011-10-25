checkSignals<-function(CNOlist, Model){

#check that CNOlist is a CNOlist
	if(!is.list(CNOlist)){
		stop("This function expects as input a CNOlist as output by makeCNOlist or normaliseCNOlist")
		}
	if(all(names(CNOlist) != 
			c("namesCues",
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
		
#check that Model is a Model list
	if(!is.list(Model)){
		stop("This function expects as input a Model as output by readSif")
		}
	if(all(names(Model) != c("reacID", "namesSpecies","interMat","notMat"))){
		stop("This function expects as input a Model as output by readSif")
		}	
		
#check that all of the signals in CNOlist$namesSignals match to one species in Model$namesSpecies
	signalsMatch<-match(CNOlist$namesSignals,Model$namesSpecies,nomatch=0)
	if(any(signalsMatch == 0)){
		warning(paste
			(
			"The following signals from your CNOlist do not match any of the species in your model, and should be removed:",
			toString(CNOlist$namesSignals[which(signalsMatch == 0)])
			)
			)
		}
		
	}


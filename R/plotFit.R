plotFit<-function(OptRes){
	par(mfrow=c(2,1),mar=c(0.5,4,4,0))
	plot(
		x=OptRes$Results[,"Generation"],
		y=OptRes$Results[,"Avg_Score_Gen"],
		xlab=NA,
		xaxt="n",
		ylab="Average score of generation",
		type="l")
	par(mar=c(4,4,0,0))
	plot(
		x=OptRes$Results[,"Generation"],
		y=OptRes$Results[,"Best_score_Gen"],
		xlab="Generations",
		ylab="Best Score",
		type="l")
	}


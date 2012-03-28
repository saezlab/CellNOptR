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
#  CNO website: http://www.ebi.ac.uk/saezrodriguez/cno
#
##############################################################################
# $Id: plotCNOlistLarge.R 595 2012-02-22 17:21:47Z cokelaer $
#This function is a variant of plotCNOlist that works for bigger datasets and allows one 
#to split the plots into n=nsplit plots (across the conditions dimension)

plotCNOlistLarge<-function(CNOlist,nsplit=4, newDevice=FALSE){

#check that CNOlist is a CNOlist
	if(!is.list(CNOlist)){
		stop("This function expects as input a CNOlist as output by makeCNOlist or normaliseCNOlist")
		}
	if(
		all(names(CNOlist) != c(
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

    L = dim(CNOlist$valueCues)[1] # number of signals/Cues
    if (nsplit<=0 | nsplit>=L){
        stop("nsplit must be strictly positive and smaller than the number of signals")
    }


	splits<-dim(CNOlist$valueCues)[1]/nsplit
	splits<-floor(splits)
	CNOlistOriginal<-CNOlist
	
	for(i in 1:nsplit){
        if (newDevice==TRUE) { # this flag allows to plot the results in nsplit different 
                               # devices, which is useful for scripting, otherwise only  
                               # the last plot is shown
            dev.new()
        }
		CNOlist<-CNOlistOriginal
		
		if(i == nsplit){
			indices<-(indices[length(indices)]+1):dim(CNOlistOriginal$valueCues)[1]
			}
			
		indices<-((1:splits)+((i-1)*splits))
		CNOlist$valueCues<-CNOlist$valueCues[indices,]
		CNOlist$valueStimuli<-CNOlist$valueStimuli[indices,]
		CNOlist$valueInhibitors<-CNOlist$valueInhibitors[indices,]
		CNOlist$valueSignals$t0<-CNOlist$valueSignals$t0[indices,]
		
		for(n in 2:length(CNOlist$valueSignals)){
			CNOlist$valueSignals[[n]]<-CNOlist$valueSignals[[n]][indices,]
			}
			
		par(
			mfrow=c(nr=dim(CNOlist$valueSignals[[1]])[1]+1,nc=dim(CNOlist$valueSignals[[1]])[2]+2),
			cex=0.5,
			pch=20,
			mar=c(0.5,0.5,0,0),
			oma=c(3,2,2,2))
		yMax<-max(unlist(lapply(CNOlist$valueSignals,function(x) max(x, na.rm=TRUE))))
		yMin<-min(unlist(lapply(CNOlist$valueSignals,function(x) min(x, na.rm=TRUE))))
		xVal<-CNOlist$timeSignals
		
		for(c in 1:dim(CNOlist$valueSignals[[1]])[2]){
			par(fg="blue",mar=c(0.5,0.5,0.7,0))


            tryCatch(
    			plot(x=xVal, y=rep(-5,length(xVal)), ylim=c(yMin, yMax),xlab=NA,ylab=NA,xaxt="n",yaxt="n"),
                error=function(e) stop(c("CNOError: an error occurred while
plotting. You may want to try again with a larger number of figures nsplit
by providing the nsplit option (suggestion: nsplit=",floor(L/25), ")")))

			text(
				labels=CNOlist$namesSignals[c],
				x=((xVal[length(xVal)]-xVal[1])/2),
				y=(yMin+((yMax-yMin)/2)),
				cex=1)
			}
     
		plot(
			x=xVal, 
			y=rep(-5,length(xVal)), 
			ylim=c(yMin, yMax),
			xlab=NA,ylab=NA,xaxt="n",yaxt="n")
		text(
			labels="Stim",
			x=((xVal[length(xVal)]-xVal[1])/2),
			y=(yMin+((yMax-yMin)/2)),cex=1)
		plot(
			x=xVal, y=rep(-5,length(xVal)), 
			ylim=c(yMin, yMax),
			xlab=NA,ylab=NA,xaxt="n",yaxt="n")
		text(
			labels="Inh",
			x=((xVal[length(xVal)]-xVal[1])/2),
			y=(yMin+((yMax-yMin)/2)),cex=1)
		par(fg="black",mar=c(0.5,0.5,0,0))	
		
		for(r in 1:dim(CNOlist$valueSignals[[1]])[1]){
		
			for(c in 1:dim(CNOlist$valueSignals[[1]])[2]){
				yVal<-lapply(CNOlist$valueSignals,function(x) {x[r,c]})
				
				plot(x=xVal,y=yVal,ylim=c(yMin, yMax),xlab=NA,ylab=NA,xaxt="n",yaxt="n")
				lines(x=xVal,y=yVal,ylim=c(yMin, yMax),xlab=NA,ylab=NA,xaxt="n",yaxt="n")
				
				#add the axis annotations: if we're on the last row, add the x axis
				if(r == dim(CNOlist$valueSignals[[1]])[1]){
					axis(side=1,at=xVal)
					}	
					
				#add the axis annotations: if we're on the first column, add the y axis	
				if(c == 1){
					axis(side=2,at=c(-0.5,0,0.5))
					}
					
					
				}
				
			image(
				t(matrix(CNOlist$valueStimuli[r,],nrow=1)),
				col=c("white","black"),xaxt="n",yaxt="n")
			
			if(r == dim(CNOlist$valueSignals[[1]])[1]){
				axis(
					side=1,
					at=seq(from=0, to=1,length.out=length(CNOlist$namesStimuli)),
					labels=CNOlist$namesStimuli,las=3,cex.axis=1)
				}
			
			image(
				t(matrix(CNOlist$valueInhibitors[r,],nrow=1)),
				col=c("white","black"),xaxt="n",yaxt="n")	
			
			if(r == dim(CNOlist$valueSignals[[1]])[1]){
				axis(
					side=1,
					at=seq(from=0, to=1,length.out=length(CNOlist$namesInhibitors)),
					labels=paste(CNOlist$namesInhibitors,"-i",sep=""),
					las=3,cex.axis=1)
			}
			
			}
		}			
	}


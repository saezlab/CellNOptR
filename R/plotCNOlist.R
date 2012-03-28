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
# $Id: plotCNOlist.R 592 2012-02-22 17:18:16Z cokelaer $
plotCNOlist<-function(CNOlist){

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
		
#set graphic parameters		
	par(
		mfrow=c(nr=dim(CNOlist$valueSignals[[1]])[1]+1,nc=dim(CNOlist$valueSignals[[1]])[2]+1),
		cex=0.5,pch=20,mar=c(0.5,0.5,0,0),oma=c(3,2,2,2))
	yMax<-max(unlist(lapply(CNOlist$valueSignals,function(x) max(x, na.rm=TRUE))))
	yMin<-min(unlist(lapply(CNOlist$valueSignals,function(x) min(x, na.rm=TRUE))))
	xVal<-CNOlist$timeSignals
	
	for(c in 1:dim(CNOlist$valueSignals[[1]])[2]){
		par(fg="blue",mar=c(0.5,0.5,0.7,0))
		plot(
			x=xVal, 
			y=rep(-5,length(xVal)), 
			ylim=c(yMin, yMax),xlab=NA,ylab=NA,xaxt="n",yaxt="n")
		text(
			labels=CNOlist$namesSignals[c],
			x=((xVal[length(xVal)]-xVal[1])/2),
			y=(yMin+((yMax-yMin)/2)),cex=2)
		}
		
		plot(
			x=xVal, 
			y=rep(-5,length(xVal)), 
			ylim=c(yMin, yMax),xlab=NA,ylab=NA,xaxt="n",yaxt="n")
		text(
			labels="Cues",
			x=((xVal[length(xVal)]-xVal[1])/2),
			y=(yMin+((yMax-yMin)/2)),cex=2)
			
	par(fg="black",mar=c(0.5,0.5,0,0))
	
	for(r in 1:dim(CNOlist$valueSignals[[1]])[1]){
	
		for(c in 1:dim(CNOlist$valueSignals[[1]])[2]){
			yVal<-lapply(CNOlist$valueSignals,function(x) {x[r,c]})
			
			plot(x=xVal,y=yVal,ylim=c(yMin, yMax),xlab=NA,ylab=NA,xaxt="n",yaxt="n")
			lines(x=xVal,y=yVal,ylim=c(yMin, yMax),xlab=NA,ylab=NA,xaxt="n",yaxt="n")
			
			#add the annotation of the axis: if we're on the last row we need an x-axis
			if(r == dim(CNOlist$valueSignals[[1]])[1]){
				axis(side=1,at=xVal)
				}	
			
			#add the annotation of the axis: if we're on the first column we need a y-axis
			if(c == 1){
				axis(side=2,at=c(-0.5,0,0.5))
				}
				
			}

            # The image (cues) last columns
			if (length(CNOlist$namesInhibitors) != 0){
                data = t(matrix(c(CNOlist$valueStimuli[r,],CNOlist$valueInhibitors[r,]),nrow=1))
                # create the color vector
                if (all(data==1)==TRUE){
                    col=c("black")
                }
                else if (all(data==0)==TRUE){
                    col=c("white")
                }
                else{
                    col=c("white", "black")
                }
                image(data,col=col,xaxt="n",yaxt="n")
            }
            else{ # special case of no inhibitors
                data = t(matrix(CNOlist$valueStimuli[r,],nrow=1))
                # create the color vector
                if (all(data==1)==TRUE){
                    col=c("black")
                }
                else if (all(data==0)==TRUE){
                    col=c("white")
                }
                else{
                    col=c("white", "black")
                }
                image(data, col=col,xaxt="n",yaxt="n")
            }

        # the axis last column
        if(r == dim(CNOlist$valueSignals[[1]])[1]){
            # special case of no inhibitors
            if (length(CNOlist$namesInhibitors) == 0){
                labels = c(CNOlist$namesStimuli)
            }
            else{
                labels=c(CNOlist$namesStimuli,paste(CNOlist$namesInhibitors,"-i",sep=""))
            }

            axis(
                side=1,
                at=seq(from=0, to=1,length.out=length(CNOlist$namesCues)),
                labels=labels, las=3,cex.axis=1)   
        }
    }
}


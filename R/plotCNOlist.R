#
#  This file is part of the CNO software
#
#  Copyright (c) 2011-2012 - EMBL - European Bioinformatics Institute
#
#  File author(s): CNO developers (cno-dev@ebi.ac.uk)
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  CNO website: http://www.cellnopt.org
#
##############################################################################
# $Id: plotCNOlist.R 3635 2013-05-28 09:21:46Z cokelaer $
plotCNOlist<-function(CNOlist){

#check that CNOlist is a CNOlist

    if ((class(CNOlist)=="CNOlist")==FALSE){
         CNOlist = CellNOptR::CNOlist(CNOlist)
    }

    # save current par() and recall at end to avoid permanent changes
    oldPar = par(no.readonly = TRUE)

    #set graphic parameters
    par(
        mfrow=c(nr=dim(CNOlist@signals[[1]])[1]+1,nc=dim(CNOlist@signals[[1]])[2]+1),
        cex=0.5,pch=20,mar=c(0.5,0.5,0,0),oma=c(3,2,2,2)
    )

    yMax<-max(unlist(lapply(CNOlist@signals,function(x) max(x, na.rm=TRUE))))
    if (yMax<=1){
        yMax = 1
    }
    yMin<-min(unlist(lapply(CNOlist@signals,function(x) min(x, na.rm=TRUE))))
    xVal<-CNOlist@timepoints

    for(c in 1:dim(CNOlist@signals[[1]])[2]){
        par(fg="blue",mar=c(0.5,0.5,0.7,0))
        plot(
            x=xVal,
            y=rep(-5,length(xVal)),
            ylim=c(yMin, yMax),xlab=NA,ylab=NA,xaxt="n",yaxt="n")
        text(
            labels=colnames(CNOlist@signals[[1]])[c],
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

    for(r in 1:dim(CNOlist@signals[[1]])[1]){

        for(c in 1:dim(CNOlist@signals[[1]])[2]){
            yVal<-lapply(CNOlist@signals,function(x) {x[r,c]})
            plot(x=xVal,y=yVal,ylim=c(yMin, yMax),xlab=NA,ylab=NA,xaxt="n",yaxt="n")
            lines(x=xVal,y=yVal,ylim=c(yMin, yMax),xlab=NA,ylab=NA,xaxt="n",yaxt="n")
            tryCatch(
                {ystd = unlist(lapply(CNOlist@variances,function(x) {x[r,c]}))
                .error.bar(unlist(xVal), unlist(yVal), ystd)},
                error=function(e){}
            )

            #add the annotation of the axis: if we're on the last row we need an x-axis
            if(r == dim(CNOlist@signals[[1]])[1]){
                axis(side=1,at=xVal)
                }

            #add the annotation of the axis: if we're on the first column we need a y-axis
            if(c == 1){
                axis(side=2,at=c(-0.5,0,0.5))
                }

            }

            # The image (cues) last columns
            if (length(CNOlist@inhibitors) != 0){
                data = t(matrix(c(CNOlist@stimuli[r,],CNOlist@inhibitors[r,]),nrow=1))
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
                data = t(matrix(CNOlist@stimuli[r,],nrow=1))
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
        if(r == dim(CNOlist@signals[[1]])[1]){
            # special case of no inhibitors
            if (length(colnames(CNOlist@inhibitors)) == 0){
                labels = colnames(CNOlist@stimuli)
            }
            else{
                labels=c(colnames(CNOlist@stimuli),paste(colnames(CNOlist@inhibitors),"-i",sep=""))
            }

            axis(
                side=1,
                at=seq(from=0, to=1,length.out=length(colnames(CNOlist@cues))),
                labels=labels, las=3,cex.axis=1)
        }
    }
    par(oldPar)

}


.error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
     if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
         stop("vectors must be same length")

     positives = upper > 0
     x = x[positives]
     y = y[positives]
     upper = upper[positives]
     lower = lower[positives]
     arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
     }


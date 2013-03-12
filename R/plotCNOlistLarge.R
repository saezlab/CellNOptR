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
#  CNO website: http://www.ebi.ac.uk/saezrodriguez/cno
#
##############################################################################
# $Id: plotCNOlistLarge.R 3155 2013-01-09 15:24:58Z cokelaer $
#This function is a variant of plotCNOlist that works for bigger datasets and allows one
#to split the plots into n=nsplit plots (across the conditions dimension)

plotCNOlistLarge<-function(CNOlist,nsplit=4, newDevice=FALSE){

    if ((class(CNOlist)=="CNOlist")==FALSE){
         CNOlist = CellNOptR::CNOlist(CNOlist)
     }

    opar = par(no.readonly = TRUE)
    on.exit(par(opar))


    L = dim(CNOlist@cues)[1] # number of signals/Cues
    if (nsplit<=0 | nsplit>=L){
        stop("nsplit must be strictly positive and smaller than the number of signals")
    }


    splits<-dim(CNOlist@cues)[1]/nsplit
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
            indices<-(indices[length(indices)]+1):dim(CNOlistOriginal@cues)[1]
            }

        indices<-((1:splits)+((i-1)*splits))
        CNOlist@cues<-CNOlist@cues[indices,]
        CNOlist@stimuli<-CNOlist@stimuli[indices,]
        CNOlist@inhibitors<-CNOlist@inhibitors[indices,]
        CNOlist@signals[[1]]<-CNOlist@signals[[1]][indices,]

        for(n in 2:length(CNOlist@signals)){
            CNOlist@signals[[n]]<-CNOlist@signals[[n]][indices,]
            }

        par(
            mfrow=c(nr=dim(CNOlist@signals[[1]])[1]+1,nc=dim(CNOlist@signals[[1]])[2]+2),
            cex=0.5,
            pch=20,
            mar=c(0.5,0.5,0,0),
            oma=c(3,2,2,2))
        yMax<-max(unlist(lapply(CNOlist@signals,function(x) max(x, na.rm=TRUE))))
        yMin<-min(unlist(lapply(CNOlist@signals,function(x) min(x, na.rm=TRUE))))
        xVal<-CNOlist@timepoints

        for(c in 1:dim(CNOlist@signals[[1]])[2]){
            par(fg="blue",mar=c(0.5,0.5,0.7,0))


            tryCatch(
                plot(x=xVal, y=rep(-5,length(xVal)), ylim=c(yMin, yMax),xlab=NA,ylab=NA,xaxt="n",yaxt="n"),
                error=function(e) stop(c("CNOError: an error occurred while
plotting. You may want to try again with a larger number of figures nsplit
by providing the nsplit option (suggestion: nsplit=",floor(L/25), ")")))

            text(
                labels=colnames(CNOlist@signals[[1]])[c],
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

        for(r in 1:dim(CNOlist@signals[[1]])[1]){

            for(c in 1:dim(CNOlist@signals[[1]])[2]){
                yVal<-lapply(CNOlist@signals,function(x) {x[r,c]})

                plot(x=xVal,y=yVal,ylim=c(yMin, yMax),xlab=NA,ylab=NA,xaxt="n",yaxt="n")
                lines(x=xVal,y=yVal,ylim=c(yMin, yMax),xlab=NA,ylab=NA,xaxt="n",yaxt="n")

                #add the axis annotations: if we're on the last row, add the x axis
                if(r == dim(CNOlist@signals[[1]])[1]){
                    axis(side=1,at=xVal)
                    }

                #add the axis annotations: if we're on the first column, add the y axis
                if(c == 1){
                    axis(side=2,at=c(-0.5,0,0.5))
                    }


                }

            image(
                t(matrix(CNOlist@stimuli[r,],nrow=1)),
                col=c("white","black"),xaxt="n",yaxt="n")

            if(r == dim(CNOlist@signals[[1]])[1]){
                axis(
                    side=1,
                    at=seq(from=0, to=1,length.out=length(CNOlist@stimuli)),
                    labels=CNOlist@stimuli,las=3,cex.axis=1)
                }

            image(
                t(matrix(CNOlist@inhibitors[r,],nrow=1)),
                col=c("white","black"),xaxt="n",yaxt="n")

            if(r == dim(CNOlist@signals[[1]])[1]){
                axis(
                    side=1,
                    at=seq(from=0, to=1,length.out=length(CNOlist@inhibitors)),
                    labels=paste(CNOlist@inhibitors,"-i",sep=""),
                    las=3,cex.axis=1)
            }

            }
        }
    }


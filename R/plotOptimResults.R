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
# $Id: plotOptimResults.R 2247 2012-08-28 16:52:04Z cokelaer $
plotOptimResults<-function(
    simResults=simResults,
    expResults=expResults,
    times=times,
    namesCues=namesCues,
    namesSignals=namesSignals,
    valueCues=valueCues, formalism="new"){

    oldPar = par(no.readonly = TRUE)
    on.exit(par(oldPar))

#Set graphical parameters
    par(
        mfrow=c(nr=dim(simResults[[1]])[1]+1,nc=dim(simResults[[1]])[2]+1),
        cex=0.5,
        pch=20,
        mar=c(0.5,0.5,0,0),
        oma=c(3,2,2,2))
    yMax<-max(max(unlist(lapply(expResults,function(x) max(x, na.rm=TRUE)))),1)
    yMin<-min(min(unlist(lapply(expResults,function(x) min(x, na.rm=TRUE)))),0)
    xVal<-times

#Write the boxes on top with the labels of the signals

    for(c in 1:dim(expResults[[1]])[2]){
        par(fg="blue",mar=c(0.5,0.5,0.5,0))
        plot(
            x=xVal,
            y=rep(-5,length(xVal)),
            ylim=c(yMin, yMax),
            xlab=NA,
            ylab=NA,
            xaxt="n",
            yaxt="n")
        text(
            labels=as.character(namesSignals[c]),
            x=((xVal[length(xVal)]-xVal[1])/2),
            y=(yMin+((yMax-yMin)/2)),
            cex=2)
        }

        plot(
            x=xVal,
            y=rep(-5,length(xVal)),
            ylim=c(yMin, yMax),
            xlab=NA,
            ylab=NA,
            xaxt="n",
            yaxt="n")

        text(
            labels="Cues",
            x=((xVal[length(xVal)]-xVal[1])/2),
            y=(yMin+((yMax-yMin)/2)),
            cex=2)

    par(fg="black",mar=c(0.5,0.5,0,0))
#Go through each elements of the results matrix
#(i.e. one plot for each signal (column) for each condition (row)

    for(r in 1:dim(expResults[[1]])[1]){

        for(c in 1:dim(expResults[[1]])[2]){

#Determine the simulated and experimental values
            yVal<-lapply(expResults,function(x) {x[r,c]})
            yValS<-lapply(simResults,function(x) {x[r,c]})

#Take care of the NA's which would mess up our mean difference used for bg colouring

            if(sum(is.na(yVal)) != 0){
                yVal4Diff<-yVal[-which(is.na(yVal))]
                yValS4Diff<-yValS[-which(is.na(yVal))]
                }else{
                    yVal4Diff<-yVal
                    yValS4Diff<-yValS
                    }

            if(sum(is.na(yValS4Diff)) != 0){
                yVal4Diff<-yVal4Diff[-which(is.na(yValS4Diff))]
                yValS4Diff<-yValS4Diff[-which(is.na(yValS4Diff))]
                }else{
                    yVal4Diff<-yVal4Diff
                    yValS4Diff<-yValS4Diff
                    }

#Compute the mean difference between data and simulation, taking into account t0 but not NAs
            if (formalism == "new"){
                diff<-mean(abs(unlist(yVal4Diff)[1:length(yVal4Diff)]-unlist(yValS4Diff)[1:length(yValS4Diff)]),na.rm=TRUE)
            }
            else{
                diff<-mean(abs(unlist(yVal4Diff)[2:length(yVal4Diff)]-unlist(yValS4Diff)[2:length(yValS4Diff)]),na.rm=TRUE)
           } 

#Set the bg colour based on the above
            if(is.na(diff)){
                bgcolor<-'white'
                }else{
                    if(diff > 0.9){
                        bgcolor<-'red'
                        }else{
                            if(diff > 0.8){
                                bgcolor<-'indianred1'
                                }else{
                                    if(diff > 0.7){
                                        bgcolor<-'lightpink2'
                                        }else{
                                            if(diff > 0.6){
                                                bgcolor<-'lightpink'
                                                }else{
                                                    if(diff > 0.5){
                                                        bgcolor<-'mistyrose'
                                                        }else{
                                                            if(diff > 0.4){
                                                                bgcolor<-'palegoldenrod'
                                                                }else{
                                                                    if(diff > 0.3){
                                                                        bgcolor<-'palegreen'
                                                                        }else{
                                                                            if(diff > 0.2){
                                                                                bgcolor<-'darkolivegreen3'
                                                                                }else{
                                                                                    if(diff > 0.1){
                                                                                        bgcolor<-'chartreuse3'
                                                                                        }else{
                                                                                            bgcolor<-'forestgreen'
                                                                                            }
                                                                                    }
                                                                            }
                                                                    }
                                                            }
                                                    }
                                            }
                                    }
                            }
                    }

#Plot
            plot(x=xVal,y=yVal,ylim=c(yMin, yMax),xlab=NA,ylab=NA,xaxt="n",yaxt="n")
            rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = bgcolor)
            points(x=xVal,y=yVal,ylim=c(yMin, yMax),xlab=NA,
                    ylab=NA,xaxt="n",yaxt="n")
            lines(x=xVal,y=yVal,ylim=c(yMin, yMax),xlab=NA,
                    ylab=NA,xaxt="n",yaxt="n")
            lines(x=xVal,y=yValS,ylim=c(yMin, yMax),xlab=NA,
                    ylab=NA,xaxt="n",yaxt="n",col="blue",lty=2)
            points(x=xVal,y=yValS,ylim=c(yMin, yMax),xlab=NA,
                    ylab=NA,xaxt="n",yaxt="n",col="blue")

#add the axis annotations: if we're on the last row, add the x axis
            if(r == dim(expResults[[1]])[1]){
                axis(side=1,at=xVal)
                }

#add the axis annotations: if we're on the first column, add the y axis
            if(c == 1){
                axis(side=2,at=c(-0.5,0,0.5))
                }


            }
    #Plot the image plots that tell preseence/absence of cues
    #
    data = t(matrix(as.numeric(valueCues[r,]),nrow=1))
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

        if(r == dim(expResults[[1]])[1]){
            axis(
                side=1,
                at=seq(from=0, to=1,length.out=length(namesCues)),
                labels=as.character(namesCues),
                las=3)
            }
        }
}


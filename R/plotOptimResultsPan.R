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
# $Id: plotOptimResultsPan.R 802 2012-03-22 16:44:12Z cokelaer $

plotOptimResultsPan <- function(simResults, yInterpol=NULL, xCoords=NULL,
CNOlist=CNOlist, formalism=c("ss1","ss2","ssN","dt","ode"), pdf=FALSE,
pdfFileName="", tPt=NULL, plotParams=list(margin=0.1, width=15, height=12,
cmap_scale=1, cex=1.6, ymin=NULL, F=1, rotation=0)) {


    if ((class(CNOlist)=="CNOlist")==FALSE){
         CNOlist = CellNOptR::CNOlist(CNOlist)
     }

    formalism <- match.arg(formalism)
    # index valueSignals according to tPt
    valueSignalsI = sapply(c(0,tPt), function(x) which(CNOlist@timepoints==x))

    if ("margin" %in% names(plotParams) == FALSE){
        plotParams$margin = 0.1
    }
    if ("cex" %in% names(plotParams) == FALSE){
        plotParams$cex = 1.6
    }
    if ("width" %in% names(plotParams) == FALSE){
        plotParams$width = 15
    }
    if ("height" %in% names(plotParams) == FALSE){
        plotParams$height = 12
    }
    if ("cmap_scale" %in% names(plotParams) == FALSE){
        plotParams$cmap_scale = 1
    }
    if ("lwd" %in% names(plotParams) == FALSE){
        plotParams$lwd = 3
    }
    if ("pch" %in% names(plotParams) == FALSE){
        plotParams$pch = 16
    }
    if ("ymin" %in% names(plotParams) == FALSE){
        plotParams$ymin = NULL
    }
    if ("F" %in% names(plotParams) == FALSE){
        plotParams$F = 1
    }
    if ("rotation" %in% names(plotParams) == FALSE){
        plotParams$rotation = 0
    }
    
# aliases
    margin = plotParams$margin
    cex = plotParams$cex

    #####    functions    #####

    list2Array = function(x, dim) {
        xUnlist = unlist(x);
        xArray=array(xUnlist,dim=dim)
    }

    #####    /functions/    #####


    opar = par(no.readonly = TRUE)
    on.exit(par(opar))

    if(pdf==TRUE) {
        pdf(file=pdfFileName, width=plotParams$width, height=plotParams$height)
    }

    Ncols = dim(CNOlist@signals[[1]])[2]
    Nrows = dim(CNOlist@signals[[1]])[1]

    split.screen(c(1, Ncols+3))
    for(a in 1:(Ncols+2)) {
        split.screen(c(Nrows+1,1),a) # why +2 ?
    }


    # TODO - do i need all these with split.screen?
    par(
        pch=2,
        # margin bottom of pcolor mesh: oma(X,Y,W,Z) where X is distance bottom
        # to the colormesh
        oma=c(3,2,2,2),
        mgp=c(3,0.9,0),
        family="Times"
    )

    Ncolors = 1000

    heatCols = heat.colors(Ncolors)

    # maximum across all data points
    #yMax <- max(unlist(lapply(CNOlist$valueSignals, function(x) max(x,na.rm=TRUE))))
    yMax=1
    # minimum across all data points

    if (is.null(plotParams$ymin)==FALSE){
        yMin = plotParams$ymin
    } else {
        yMin <- min(unlist(lapply(CNOlist@signals, function(x) min(x,na.rm=TRUE))))
    }


    # time labels
    xVal <- CNOlist@timepoints[valueSignalsI]
    if(formalism=="dt") {
        xValS = xCoords
        norm = length(CNOlist@timepoints)
    } else if (formalism == "ss1") {
        xValS = c(0,tPt[1])
        norm = 2
    } else if (formalism == "ss2") {
        xValS = c(0,tPt[1:2])
        norm = 3
    } else if (formalism == "ssN") {
      xValS = c(0,tPt)
      norm = length(tPt)+1 # should be number of time points.
    } else if (formalism == "ode") {
        xValS = CNOlist@timepoints
        xVal <- CNOlist@timepoints
        valueSignalsI = seq(1, length(xValS))
        norm = length(CNOlist@timepoints)
    }else {
        xValS = xVal
    }
    # latest time point
    xValMax = max(xVal)

    # make simResults array if not already

    if(!is.array(simResults)) {
        simResults = list2Array(simResults, dim=c(dim(simResults[[1]]),length(simResults)))
    }

    # make valueSignals an array
    valueSignalsArr = list2Array(CNOlist@signals,
        dim=c(dim(CNOlist@signals[[1]]),length(CNOlist@signals)))

    # calculate the MSE
    allDiff = matrix(NA, nrow=dim(simResults)[1], ncol=dim(simResults)[2])


    if(formalism != "dt") {
        # ss1, ss2, ssN
        for(a in 1:dim(simResults)[1]) {
            for(b in 1:dim(simResults)[2]) {
                c1 = simResults[a,b,]
                c2 = valueSignalsArr[a,b,valueSignalsI]
                NAcount = max(sum(is.na(c1)), sum(is.na(c2)))
                allDiff[a,b] = sum((c1 - c2)^2, na.rm=T)/(norm - NAcount)
            }
        }
    } else {
        # dt 
        for(a in 1:dim(simResults)[1]) {
            for(b in 1:dim(simResults)[2]) {
                c1 =  simResults[a,b,]
                c2 = yInterpol[a,b,]
                NAcount = max(sum(is.na(c1)), sum(is.na(c2)))
                allDiff[a,b] = sum((c1 - c2)^2, na.rm=T)/(norm - NAcount)
            }
        }
    }
    # max difference between sim and exper
    #diffMax = max(unlist(!is.na(allDiff)))

    # set the count for the split screen window
    count1 = dim(CNOlist@signals[[1]])[2]+4

    # plot headers
    for(c in 1:dim(CNOlist@signals[[1]])[2]) {
        screen(count1)
        par(fg="blue",mar=c(margin, margin, margin, 0))
        plot(x=xVal, y=rep(-5,length(xVal)), ylim=c(yMin, yMax),
        xlab=NA,ylab=NA,xaxt="n",yaxt="n")
		labels=colnames(CNOlist@signals[[1]])[c]
        text(
            labels=labels,
            x=((xVal[length(xVal)]-xVal[1])/2),
            y=(yMin+((yMax-yMin)/2)),
            cex=cex
        )
        count1 = count1 + dim(CNOlist@signals[[1]])[1]+1
    }

    # stim + inhib
    screen(count1)
    par(fg="blue",mar=c(margin, margin, margin, 0))
    plot(
        x = xVal,
        y = rep(-5,length(xVal)),
        ylim = c(yMin, yMax),
        xlab = NA,ylab=NA,xaxt="n",yaxt="n"
    )
    text(
        labels = "Stim",
        x = ((xVal[length(xVal)]-xVal[1])/2),
        y = (yMin+((yMax-yMin)/2)),cex=cex
    )

    count1 = count1 + dim(CNOlist@signals[[1]])[1]+1
       screen(count1)
    #par(fg="blue",mar=c(0.5,0.5,0.7,0))
    par(fg="blue",mar=c(margin, margin, margin, 0))
       plot(
        x = xVal, y=rep(-5,length(xVal)),
        ylim = c(yMin, yMax),
           xlab = NA,ylab=NA,xaxt="n",yaxt="n")
    text(
        labels="Inh",
         x=((xVal[length(xVal)]-xVal[1])/2),
           y=(yMin+((yMax-yMin)/2)),cex=cex
       )


    # new count for plotting results
    countRow = dim(CNOlist@signals[[1]])[2]+4


    for(c in 1:dim(CNOlist@signals[[1]])[2]) {
        countRow=countRow+1
        for(r in 1:dim(CNOlist@signals[[1]])[1]) {

            screen(countRow)
            par(fg="black",mar=c(margin, margin,0,0))
            yVal <- lapply(CNOlist@signals[valueSignalsI], function(x) {x[r,c]})

            yValS <- simResults[r,c,]
            if(!is.na(allDiff[r,c])) {
                #diff = (1 - (allDiff[r,c] / diffMax)) * 1000
                diff = (1 - (allDiff[r,c] /plotParams$F)**plotParams$cmap_scale) * (Ncolors-1)+1
            } else {
                diff = -1
            }
            if (diff <1 & diff!=-1){diff=1}

            if(diff>Ncolors) {
                diff=Ncolors
            }
            # this is a NA or NaN so what color ? 
            if(diff==-1) {
                bgcolor = "gray"
            } else{
                bgcolor = heatCols[diff]
            }
            plot(x=xVal,y=yVal,ylim=c(yMin,yMax),xlab=NA,ylab=NA,xaxt="n",yaxt="n",)
            rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col=bgcolor)

            
            tryCatch(
                {ystd = unlist(lapply(CNOlist@variances[valueSignalsI],function(x) {x[r,c]}))
                .error.bar(unlist(xVal), unlist(yVal), ystd)
                }, error=function(e){}
            )




            # the measurements
            lines(xVal, yVal, ylim=c(yMin, yMax), xlab=NA, ylab=NA, xaxt="n", yaxt="n", col="black", lty=1, lwd=1)
            points(x=xVal,y=yVal,ylim=c(yMin, yMax),xlab=NA,ylab=NA,xaxt="n",yaxt="n",pch=2)

            # the simulated data
            lines(xValS, yValS, xlim=c(0,xValMax), ylim=c(yMin, yMax),
                xlab=NA, ylab=NA, xaxt="n", yaxt="n", col="blue", lty=2, lwd=plotParams$lwd)
            points(xValS, yValS, xlim=c(0,xValMax), ylim=c(yMin, yMax),
                xlab=NA, ylab=NA, xaxt="n", yaxt="n", col="blue", pch=plotParams$pch)

            # add the axis annotations: if we're on the last row, add the x axis
            if(r == dim(CNOlist@signals[[1]])[1]) {
                axis(side=1,at=CNOlist@timepoints)
            }

            # add the axis annotations: if we're on the first column, add the y axis
            if(c == 1)    {
                axis(side=2,at=c(0,0.5,1), labels=c("","0.5",""), las=1)
            }

            countRow=countRow+1
        }
    }

    sStim = countRow+1

    for(c1 in 1:dim(CNOlist@signals[[1]])[1]) {
    	screen(sStim)
    	#par(mar=c(0.5,0.5,0.5,0.5))
    	par(mar=c(margin, margin,0,0))
    	
    	
    	data = matrix(CNOlist@stimuli[c1,],nrow=1)
    	if(c1 == dim(CNOlist@signals[[1]])[1]){
    		barplot(data,yaxt="n",ylim=c(0,1),names.arg = colnames(CNOlist@stimuli),las=2)
    		#axis(1)
    	}else{
    		barplot(data,xaxt="n",yaxt="n",ylim=c(0,1))
    		#axis(1)
    	}
    	
    	
    	sStim = sStim+1
    }

    sInhib = sStim+1
    for(c1 in 1:dim(CNOlist@signals[[1]])[1]) {
        screen(sInhib)
    	
    	#par(mar=c(0.5,0.5,0.5,0.5))
    	par(mar=c(margin, margin,0,0))
    	
    	if (length(CNOlist@inhibitors) != 0){
    		data = matrix(c(CNOlist@inhibitors[c1,]),nrow=1)
    		if(c1 == dim(CNOlist@signals[[1]])[1]){
    			barplot(data,yaxt="n",ylim=c(0,1),names.arg = c(paste(colnames(CNOlist@inhibitors),"-i",sep="")),las=3 )
    			#axis(1)
    		}else{
    			barplot(data,xaxt="n",yaxt="n",ylim=c(0,1))
    			#axis(1)
    		}
    	}
    	else{ # special case of no inhibitors
    		image(
    			t(matrix(0,nrow=1)),
    			col=c("grey"),xaxt="n",yaxt="n"
    		)
    	}
    	
    	
        sInhib = sInhib+1
    }

    screen(dim(CNOlist@signals[[1]])[2]+3)
    splitProp = 1/(dim(CNOlist@signals[[1]])[1]+1)
    sSplit = matrix(c(0,1,(1-splitProp),1,0,1,0,(1-splitProp)),ncol=4, byrow=T)
    split.screen(sSplit)
    screen(sInhib)
#       par(fg="blue",mar=c(0.5,0.5,0.7,0))
    par(fg="blue",mar=c(margin, margin, margin, 0))
    plot(
        x = xVal,
        y = rep(-5,length(xVal)),
        ylim = c(yMin, yMax),
        xlab = NA,ylab=NA,xaxt="n",yaxt="n"
    )
       text(
        labels = "Error",
        x = ((xVal[length(xVal)]-xVal[1])/2),
        y = (yMin+((yMax-yMin)/2)),cex=cex
        )

    screen(sInhib+1)

    Ncolors = 100
    colbar = heat.colors(Ncolors)

    # scale the colorbar in the same way as the colors used in the boxes
    # according to diff.
    colbar = colbar[as.integer((seq(0, Ncolors,length.out=100)/Ncolors)**plotParams$cmap_scale*Ncolors)]

    labels = c(1,0.5,0)
    len <- length(colbar)
    rhs <- 0.5
    rhs2 <- rhs + rhs/10
    at <- c(0, 0.5, 1)

    par(mai=c(0,0.1,0,0))
    #plot.new()

    tryCatch({
        plot.new()
        yyy <- seq(0,1,length=len+1)
        rect(0, yyy[1:len], rep(rhs, len), yyy[-1], col = colbar, border = colbar)
        rect(0, 0, rhs, 1, border="black")
        segments(rhs, at, rhs2, at)
        text(x=rhs2, y=at, labels=labels, pos=4, offset=0.2)

        },
        error=function(e){print("could not add colorbar. Skipped")}
    )


    close.screen(all.screens=TRUE)
#    par(oldPar)
    if(pdf==TRUE) {
        dev.off()
    }

    return(allDiff)
}


.error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
     if(length(x) != length(y) | length(y) !=length(lower) | length(lower) !=
length(upper))
         stop("vectors must be same length")

     positives = upper > 0 
     x = x[positives]
     y = y[positives]
     upper = upper[positives]
     lower = lower[positives]
     arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
     }   


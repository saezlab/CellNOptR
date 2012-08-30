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
# $Id$

## given a MIDAS file, returns a CNOlist object
## You may want to make this (and the existing readMIDAS)
## into S4 generics in order to avoid a name conflict.





## Class definition
setClass("CNOlist",
    representation(
        cues="matrix",
        inhibitors="matrix",
        stimuli="matrix",
        signals="list",
        timepoints="vector"),
    ## validity method
    validity=function(object) {
        msg <- NULL
        nrow <- nrow(cues(object))
        signalNrows <- unique(sapply(signals(object), nrow))
        if (nrow != nrow(inhibitors(object)) ||
            nrow != nrow(stimuli(object)) ||
            length(signalNrows) != 1 ||
            nrow != signalNrows)
        {
            msg <- "'nrow' differs between elements"
        }
        if (is.null(msg)) TRUE else msg
    })


## constructor
CNOlist <-function(data, subfield=FALSE, verbose=FALSE){

    # input can a filename or the old (still used) CNOlist returned by
    # makeCNOlist function. subfield and verbose used only if MIDASfile is a
    # string.
    if (is.character(data)== TRUE){
        res = internal_CNOlist_from_file(data, subfield, verbose)
    }else {
        if (is.list(data)==TRUE){
            if ("namesCues" %in% names(data) == TRUE){
                res = internal_CNOlist_from_makeCNOlist(data) 
            }else{
                stop("Not a valid list. Does not seem to be returned by CellNOptR::makeCNOlist")
            }
        }else{
        stop("Input data must be a filename or the output of CellNOptR::makeCNOlist function")
        }
    }


    new("CNOlist", cues=res$cues, inhibitors=res$inhibitors,
        stimuli=res$stimuli, signals=res$signals, timepoints=res$timepoints)
}



## accessors
cues <- function(cnoList, ...) cnoList@cues
inhibitors <- function(cnoList, ...) cnoList@inhibitors
stimuli <- function(cnoList, ...) cnoList@stimuli
signals <- function(cnoList, ...) cnoList@signals
timepoints <- function(cnoList, ...) cnoList@timepoints
#plot <- function(cnoList, ...) plotCNOlist(cnoList)


## show method
setMethod(show, "CNOlist", function(object) {
    cat("class:", class(object), "\n")
    cat("cues:", colnames(cues(object)), "\n")
    cat("inhibitors:", colnames(inhibitors(object)), "\n")
    cat("stimuli:", colnames(stimuli(object)), "\n")
    cat("timepoints:", names(signals(object)), "\n")
    cat("signals:", colnames(signals(object)[[1]]), "\n")
})

setMethod("plot", signature(x="CNOlist", y="missing"), function(x, y){
    plotCNOlist(x)
})


# used by the constructor not for export.
# open a MIDAS file (given a  filename) and create the instance of CNOlist
internal_CNOlist_from_file <- function(MIDASfile, subfield=FALSE, verbose=FALSE)
{
    x <- readMIDAS(MIDASfile, verbose=verbose)
    cnolist <- makeCNOlist(x, subfield=subfield, verbose=verbose)
    res <- internal_CNOlist_from_makeCNOlist(cnolist)
    return(res)
}


# used by the constructor not for export.
# open a MIDAS file (given a  filename) and create the instance of CNOlist
internal_CNOlist_from_makeCNOlist <- function(cnolist)
{

    myCues <- cnolist$valueCues
    colnames(myCues) <- cnolist$namesCues

    myInhibitors <- cnolist$valueInhibitors
    colnames(myInhibitors) <- cnolist$namesInhibitors

    myStimuli <- cnolist$valueStimuli
    colnames(myStimuli) <- cnolist$namesStimuli

    mySignals <- cnolist$valueSignals
    names(mySignals) <- cnolist$timeSignals
    mySignals <-
        lapply(mySignals, "colnames<-", cnolist$namesSignals)

    myTimePoints <- cnolist$timeSignals

    #CNOlist(myCues, myInhibitors, myStimuli, mySignals)
    return( list(cues=myCues, inhibitors=myInhibitors, stimuli=myStimuli,
        signals=mySignals, timepoints=myTimePoints))
}


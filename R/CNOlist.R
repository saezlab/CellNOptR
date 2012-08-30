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
#CNOlist <-function(cues, inhibitors, stimuli, signals, ...){
#    new("CNOlist", cues=cues, inhibitors=inhibitors,
#        stimuli=stimuli, signals=signals, ...)
#}
CNOlist <-function(MIDASfile, subfield=FALSE, verbose=FALSE){
    res = internal_CNOlist(MIDASfile, subfield, verbose)
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








# used by the constructor not for export.
internal_CNOlist <- function(MIDASfile, subfield=FALSE, verbose=FALSE)
{
    x <- readMIDAS(MIDASfile, verbose=verbose)
    y <- makeCNOlist(x, subfield=subfield, verbose=verbose)

    myCues <- y$valueCues
    colnames(myCues) <- y$namesCues

    myInhibitors <- y$valueInhibitors
    colnames(myInhibitors) <- y$namesInhibitors

    myStimuli <- y$valueStimuli
    colnames(myStimuli) <- y$namesStimuli

    mySignals <- y$valueSignals
    names(mySignals) <- y$timeSignals
    mySignals <-
        lapply(mySignals, "colnames<-", y$namesSignals)

    myTimePoints <- y$timeSignals

    #CNOlist(myCues, myInhibitors, myStimuli, mySignals)
    return( list(cues=myCues, inhibitors=myInhibitors, stimuli=myStimuli,
        signals=mySignals, timepoints=myTimePoints))
}


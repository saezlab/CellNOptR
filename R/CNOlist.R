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
        variances="list",
        timepoints="vector"),
    ## validity method
    validity=function(object) {
        msg <- NULL
        nrow <- nrow(getCues(object))
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
CNOlist <-function(data, verbose=FALSE){

    res = NULL
    # input can a filename or the old (still used) CNOlist returned by
    # makeCNOlist function. 
    if (is.character(data)== TRUE){
        res = internal_CNOlist_from_file(data, verbose)
    }

    if (is.list(data)==TRUE){
        if ("namesCues" %in% names(data) == TRUE){
            res = internal_CNOlist_from_makeCNOlist(data) 
        }else{
            stop("Not a valid list. Does not seem to be returned by CellNOptR::makeCNOlist")
        }
    }

    if (class(data)=="CNOlist"){
        # nothing to do already a cnolist
        res = list(
                cues=data@cues,
                inhibitors=data@inhibitors,
                stimuli=data@stimuli,
                signals=data@signals,
                variances=data@variances,
                timepoints=data@timepoints)
    }
 
    if (is.null(res)==TRUE){
        stop("Input data must be a filename or the output of CellNOptR::makeCNOlist function")
    }

    new("CNOlist", cues=res$cues, inhibitors=res$inhibitors,
        stimuli=res$stimuli, signals=res$signals, variances=res$variances, timepoints=res$timepoints)
}

setGeneric("getCues", function(object){standardGeneric("getCues")})
setGeneric("getInhibitors", function(object){standardGeneric("getInhibitors")})
setGeneric("getStimuli", function(object){standardGeneric("getStimuli")})
setGeneric("getSignals", function(object){standardGeneric("getSignals")})
setGeneric("getVariances", function(object){standardGeneric("getVariances")})
setGeneric("getTimepoints", function(object){standardGeneric("getTimepoints")})


setMethod("getCues", "CNOlist", function(object){return(object@cues)})
setMethod("getInhibitors", "CNOlist", function(object){return(object@inhibitors)})
setMethod("getStimuli", "CNOlist", function(object){return(object@stimuli)})
setMethod("getSignals", "CNOlist", function(object){return(object@signals)})
setMethod("getVariances", "CNOlist", function(object){return(object@variances)})
setMethod("getTimepoints", "CNOlist", function(object){
    # timepoints may have been modify on the fly so let us recomput it
    object@timepoints = names(object@signals)
    return(object@timepoints)

})


# timepoints will be updated if signals is changed so we should not provide any
# setTimepoint method


setGeneric("setSignals<-",function(object,value){standardGeneric("setSignals<-")})
setReplaceMethod("setSignals","CNOlist", 
    function(object,value){
        object@signals<-value
        return(object)
    }
)

## internal accessors only ??
#cues <- function(cnoList, ...) cnoList@cues
inhibitors <- function(cnoList, ...) cnoList@inhibitors
stimuli <- function(cnoList, ...) cnoList@stimuli
signals <- function(cnoList, ...) cnoList@signals


## show method
setMethod(show, "CNOlist", function(object) {
    cat("class:", class(object), "\n")
    cat("cues:", colnames(getCues(object)), "\n")
    cat("inhibitors:", colnames(getInhibitors(object)), "\n")
    cat("stimuli:", colnames(getStimuli(object)), "\n")
    cat("timepoints:", names(signals(object)), "\n")
    cat("signals:", colnames(signals(object)[[1]]), "\n")
    cat("variances:", colnames(signals(object)[[1]]), "\n")
    cat("--\nTo see the values of any data contained in this instance, just use the
appropriate getter method (e.g., getCues(cnolist), getSignals(cnolist), ...\n\n")
})

#setMethod("plot", signature(x="CNOlist", y="missing"), function(x, y, ...){
#    plotCNOlist(x)
#})

setMethod("plot", "CNOlist", function(x, y, ... ){
    plotCNOlist(x)
})
setMethod("plot", signature(x="CNOlist", y="CNOlist"), function(x, y, ... ){
    plotCNOlist2(x,y)
})

setMethod("length", "CNOlist", function(x) length(x@signals))

if (isGeneric("randomize")==FALSE){
    setGeneric(
        name="randomize",
        def=function(object,sd=0.1, minValue=0,maxValue=1,mode="gaussian"){standardGeneric("randomize")}
    )
}
#lockBinding("randomize", .GlobalEnv)


setMethod("randomize", "CNOlist", 
    definition=function(object, sd=0.1, minValue=0, maxValue=1,mode="uniform"){
        res = randomizeCNOlist(object, sd=sd, mode=mode)
        return(res)
    }
)

# a method that convert back the CNOlist to a makeCNOlist useful for back
# compatibility with other packages (e.g. ODE)
setGeneric(
    name="compatCNOlist",
    def=function(object){standardGeneric("compatCNOlist")})
setMethod("compatCNOlist", "CNOlist",
    definition=function(object){
        return(internal_compatCNOlist(object))})


internal_compatCNOlist<-function(cnolist){


    if (class(cnolist)=="CNOlist"){

        # conversion
        res = list(
             namesCues=colnames(cnolist@cues),
             namesStimuli=colnames(cnolist@stimuli),
             namesInhibitors=colnames(cnolist@inhibitors),
             namesSignals=colnames(cnolist@signals[[1]]),
             timeSignals=getTimepoints(cnolist),
             valueCues=cnolist@cues,
             valueInhibitors=cnolist@inhibitors,
             valueStimuli=cnolist@stimuli,
             valueSignals=cnolist@signals)

    } else{
        res = cnolist
    }

    return(res)
}



# used by the constructor not for export.
# open a MIDAS file (given a  filename) and create the instance of CNOlist
internal_CNOlist_from_file <- function(MIDASfile, verbose=FALSE)
{
    x <- readMIDAS(MIDASfile, verbose=verbose)

    # in the old makeCNOlist, there is a need for subfield. In cnolist, we
    # automatically figure it out here below:
    names = colnames(x$dataMatrix[x$TRcol])
    
    # Are all TR coded using subfield ?
    all_subfield_s = sapply(names, function(x){x = grepl(":Stimuli", x)})
    all_subfield_i = sapply(names, function(x){x = grepl(":Inhibitors", x)})
    all_subfield = all(all_subfield_s | all_subfield_i)

    # any subfield used ? 
    any_subfield_s = any(sapply(names, function(x){x = grepl(":Stimuli", x)}))
    any_subfield_i = any(sapply(names, function(x){x = grepl(":Inhibitors", x)}))

    

    any_subfield = any_subfield_s | any_subfield_i

    if (all_subfield == F && any_subfield == T){
        stop("If you use subfield in the MIDAS header, they should be used for all treatments adding :Stimuli or :Inhibitors after the species name (e.g., TR:<specy name>:Stimuli)")
    }

    if (all_subfield == F){
        subfield = F
    } else{
        subfield = T
    }
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
    mySignals <- lapply(mySignals, "colnames<-", cnolist$namesSignals)

    if ("valueVariances" %in% names(cnolist)){
        myVars <- cnolist$valueVariances
        myVars <- lapply(myVars, "colnames<-", cnolist$namesSignals)
    } else{
        myVars = mySignals
        for (time in 1:length(mySignals)){
            myVars[[time]] = myVars[[time]] * NA
        }
    }

    myTimePoints <- as.numeric(cnolist$timeSignals)

    #CNOlist(myCues, myInhibitors, myStimuli, mySignals)
    return( list(cues=myCues, inhibitors=myInhibitors, stimuli=myStimuli,
        signals=mySignals, variances=myVars, timepoints=myTimePoints))
}



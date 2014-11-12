readSif <- function(sifFile){
    warning("readSif is a deprecated function. Use readSIF instead. Calling readSIF for you")
    return(readSIF(sifFile))
}

prep4Sim <- function(model, params){
    warning("prep4Sim is a deprecated function. Use prep4sim instead. ")
    return(prep4sim(model))
}


simulateT1 <- function(CNOlist, model, bStringT1,simList, indexList){
    warning("simulateT1 is a deprecated function. Use simulateTN instead. ")
    return(simulateTN(CNOlist, model, bStrings=list(bStringT1)))
}

simulateT2 <- function(CNOlist, model, bStringT1,simList, indexList){
    warning("simulateT1 is a deprecated function. Use simulateTN instead. ")
    return(simulateTN(CNOlist, model, bStrings=list(bStringT1)))
}


gaBinaryT2 <- function(CNOlist, model, simList, indexList, bStringT1,
    sizeFac=0.0001, NAFac=1, popSize=50, pMutation=0.5, maxTime=60,
    maxGens=500, stallGenMax=100, selPress=1.2, elitism=5, relTol=0.1,
    verbose=TRUE, priorBitString=NULL, maxSizeHashTable=5000){

    warning("gaBinaryT2 is a deprecated function. Use gaBinaryTN instead.")
    return(gaBinaryTN(CNOlist=CNOlist, model=model, bStrings=list(bStringT1),
     sizeFac=sizeFac, NAFac=NAFac, popSize=popSize,  pMutation=pMutation,  maxTime=maxTime,
     maxGens=maxGens,  stallGenMax=stallGenMax,  selPress=selPress,  elitism=elitism,
     relTol=relTol,  verbose=verbose, priorBitString=priorBitString,
maxSizeHashTable=maxSizeHashTable))
}


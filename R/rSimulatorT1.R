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
# $Id$

# this is the code from simulatorT1 simulatorT1 was written in c

rSimulatorT1 <- function(CNOlist, model, simList, indexList) {


    if ((class(CNOlist)=="CNOlist")==FALSE){
         CNOlist = CellNOptR::CNOlist(CNOlist)
     }
    
	nSp <- dim(model$interMat)[1] 
	nReacs <- dim(model$interMat)[2]
	nCond <- dim(CNOlist@stimuli)[1]
	maxIpg <- dim(simList$finalCube)[2]

    if(is.null(dim(model$interMat))) {
		nSp <- length(model$interMat)
		nReacs <- 1
        maxIpg <- length(simList$finalCube)
	}

	# this holds, for each sp, how many reactions have that sp as output
    endIx <- rep(NA,nSp)

	for(i in 1:nSp) {
		endIx[i]<-length(which(simList$maxIx == i))
	}

	maxgpo <- max(endIx)

    # this value is used to test the stop condition for difference between 2
    # iterations.
    testVal <- 1E-3

	# create an initial values matrix
    initValues <- matrix(data=NA, nrow=nCond, ncol=nSp)
    colnames(initValues) <- model$namesSpecies

    # see warning on the top
    # set the initial values of the stimuli
    initValues[,indexList$stimulated] <- CNOlist@stimuli

    # see warning on the top
    # flip the inhibitors so that 0=inhibited/1=non-inhibitted
    valueInhibitors <- 1-CNOlist@inhibitors
    valueInhibitors[which(valueInhibitors == 1)] <- NA

    # see warning on the top
    # set the initial values of the inhibited species: 0 if inhibited, untouched
    # if not inhibited.
    initValues[,indexList$inhibited] <- valueInhibitors

    # initialise main loop
    newInput <- initValues
    termCheck1<- TRUE
    termCheck2 <- TRUE
    count <- 1

    minNA <- function(x) {
		if(all(is.na(x))) {
			return(NA)
        } else {
            return(min(x,na.rm=TRUE))
        }
    }

    filltempCube <- function(x) {
        cMatrix <- matrix(data=x,nrow=nReacs,ncol=nCond)
        cVector <- apply(cMatrix,1,function(x){return(x)})
    	return(cVector)
    }

    # the ixNeg and ignore matrices can be filled once for all before the loop
    if(nReacs > 1) {
        tempIxNeg <- apply(simList$ixNeg,2,filltempCube)
        tempIgnore <- apply(simList$ignoreCube,2,filltempCube)
    } else {
        tempIxNeg < -matrix(
        	simList$ixNeg,nrow=nCond,
            ncol=length(simList$ixNeg),byrow=TRUE)
        tempIgnore <- matrix(
            simList$ignoreCube,nrow=nCond,
            ncol=length(simList$ignoreCube),byrow=TRUE)
    }

    # main loop
    while(termCheck1 && termCheck2 ){

        outputPrev <- newInput
        
        # this is now a 2 columns matrix that has a column for each input
        # (column in finalCube)
        # and a set of rows for each reac
        # (where a set contains as many rows as conditions)
        # all concatenated into one long column.
        if(nReacs > 1) {
            tempStore <- apply(simList$finalCube,2,function(x){return(outputPrev[,x])})
        } else {
            tempStore <- outputPrev[,simList$finalCube]
        }

        # compute all the ands by taking, for each gate,
        # the min value across the inputs of that gate
        if(nReacs > 1){
            
            # new code with inf instead of NA. Speed up the code
            # set to Inf the values that are "dummies", so they won't influence the min
            tempStore[tempIgnore] <- Inf

            # flip the values that enter with a negative sign
            tempStore[tempIxNeg] <- 1-tempStore[tempIxNeg]

            # IMPORTANT NOTE: when multiple inputs and one or more is NA, then
            # this code should return NA (not the case in this frozen code). The
            # C simulator has the fix.
            outputCube <- apply(tempStore, 1, min)

            # outputCube is now a vector of length (nCond*nReacs) that contains the input of each reaction in
            # each condition, concatenated as such allcond4reac1, allcond4reac2, etc...
            # this is transformed into a matrix with a column for each reac and a row for each cond
            outputCube <- matrix(outputCube, nrow=nCond, ncol=nReacs)

            # go through each species, and if it has inputs, then take the max across those input reactions
            # i.e. compute the ORs

            # create the newInput matrix once for all. Faster than implicit
            # creation
            newInput = matrix(ncol=nSp, nrow=nCond, byrow=FALSE)
            for(s in 1:nSp){
                if(endIx[s] != 0) {
                    selection <- outputCube[,simList$maxIx == s]
                    if (length(selection) == nCond){
                        newInput[,s] <- selection
                    }
                    else{
                        selection[is.na(selection)] <- -Inf
                        newInput[,s] <- apply(selection, 1, max)
                    }
                }
            }
            newInput[is.infinite(newInput)] <- NA
        } else {
            # OLD code with NA. slower but we hardly enter in this loop
            # set to NA the values that are "dummies", so they won't influence the min
            tempStore[tempIgnore] <- NA
            # flip the values that enter with a negative sign
            tempStore[tempIxNeg] <- 1-tempStore[tempIxNeg]

            outputCube <- ifelse(all(is.na(tempStore)),NA,min(tempStore,na.rm=TRUE))
            newInput[,simList$maxIx] <- outputCube
        }

    	# reset the inhibitors and stimuli
    	for (stim in 1:length(indexList$stimulated)) {
            stimM <- cbind(
                CNOlist@stimuli[,stim],
                newInput[,indexList$stimulated[stim]])
            maxNA <- function(x){
                return(max(x,na.rm=TRUE))
            }
            stimV <- apply(stimM,1,maxNA)
            newInput[,indexList$stimulated[stim]] <- stimV
		}
    	
    	valueInhibitors <- 1-CNOlist@inhibitors
    	newInput[,indexList$inhibited] <- valueInhibitors*newInput[,indexList$inhibited]

    	# replace NAs with zeros to avoid having the NA penalty applying to unconnected species
    	newInput[is.na(newInput)] <- 0
    	outputPrev[is.na(outputPrev)] <- 0

    	# check the 2 termination conditions
    	termCheck1 <- !all(abs(outputPrev-newInput)<testVal)
   	 	termCheck2 <- (count < (nSp*1.2))
    	count <- count+1
    }

	# set the non-resolved bits to NA
    newInput[which(abs(outputPrev-newInput) > testVal)] <- NA

    return(newInput)
}


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
# $Id: cutNONC.R 4441 2014-03-10 15:01:28Z aidanmac $

cutNONC <- function(model, NONCindexes) {

    #####     FUNCTIONS    #####

    multipleInOut <- function(x) {

        spInReac = which(x != 0)
        inputs = length(which(x==-1))

        if((inputs > 1) && length(intersect(spInReac,NONCindexes))) {
            return(TRUE)
        } else {
            return(FALSE)
        }
    }

    emptyInOut<-function(x) {

        input <- match(-1,x,nomatch=0)
        output <- match(1,x,nomatch=0)

        if((input == 0) | (output == 0)) {
            return(TRUE)
        } else {
            return(FALSE)
        }
    }

    #####    /FUNCTIONS/    #####


    if(length(NONCindexes) == 0) {

        newModel=model

    } else {

        # if there are AND gates in the PKN, find out
        # if any NONCs are parts of AND gates
        editReac = apply(model$interMat,2,multipleInOut)
        newSpecies <- model$namesSpecies[-NONCindexes]
        newInterMat <- model$interMat[-NONCindexes,]
        newNotMat <- model$notMat[-NONCindexes,]

        # this function finds whether a given vector contains at least one input and one output or not
        # it returns true if the re is an in/out missing

        # modify to account for AND gates:
        # some NONC sp may be part of an AND gate with non-NONC sp
        # if this is the case, remove the NONC edge, not the reaction from reacID

        # rebuild reacIDs from matrix
        toEdit = as.numeric(which(editReac==TRUE))
        for(a in toEdit) {
            andInput = rownames(newInterMat)[which(newInterMat[,a] == -1)]
            andInd = which(newInterMat[,a] == -1)
            if(length(intersect(which(newNotMat[,a]==1),which(newInterMat[,a]==-1)))) {
                andNeg = intersect(which(newNotMat[,a]==1),which(newInterMat[,a]==-1))
                    for(p in 1:length(andNeg)) {
                        andInput[which(andInd==andNeg[p])] = paste("!", andInput[which(andInd==andNeg[p])], sep="")
                    }
            }

            LHS = paste(andInput,collapse="+", sep="")
            colnames(newInterMat)[a] = paste(LHS, "=", rownames(newInterMat)[which(newInterMat[,a] == 1)], sep="")
	    
            colnames(newNotMat)[a] = colnames(newInterMat)[a]
        }

        reac2remove <- apply(newInterMat,2,emptyInOut)
        # check if reac2remove contains selfloops
        # the condition below matches before '=' and checks this pattern after
        matchIn = gsub("^!*([[:alnum:]]+)=.*","\\1", names(reac2remove))
        potentialLoops = paste(matchIn,"=",matchIn,sep="")
        # check for + or - loops
        reac2remove[which(paste("!",potentialLoops,sep="")==names(reac2remove) | potentialLoops==names(reac2remove))] = FALSE
        
        if(any(reac2remove)) {

            reac2remove <- which(reac2remove)
            newInterMat <- newInterMat[,-reac2remove]
            newNotMat <- newNotMat[,-reac2remove]
            newreacID <- colnames(newInterMat)

            newModel <- list(
                reacID=newreacID,
                namesSpecies=newSpecies,
                interMat=newInterMat,
                notMat=newNotMat
            )
        } else {
            newModel <- list(
                reacID=colnames(newInterMat),
                namesSpecies=newSpecies,
                interMat=newInterMat,
                notMat=newNotMat)
        }
    }

    return(newModel)

}


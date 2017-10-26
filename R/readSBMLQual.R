#
#  This file is part of the CNO software
#
#  Copyright (c) 2011-2013 - EMBL-EBI
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
# $Id$

# This function reads a SBMLQual written for storing SIF file. 
# The listOfQualitativeSpecies should contain the species
# From the listOfTransitions, the edgse can be extracted
# TODO: what if there is no transittion for a given specy. Is it possible ? If
# so this is an isolated node.

readSBMLQual <- function(filename){
    warning("experimental SBML reader. use with care July 2013.")
    library(XML)
    doc = xmlTreeParse(filename)
    r = xmlRoot(doc)

    sif = data.frame(matrix(0, ncol=3))

    species = xmlChildren(xmlChildren(r)[[1]])['listOfQualitativeSpecies']
    nameSpecies = as.vector(unlist(xmlApply(species[[1]], function(x) xmlGetAttr(x, "id"))))

    transitions = xmlChildren(xmlChildren(r)[[1]])['listOfTransitions']

    andCounter = 1
    for (i in seq_along((transitions[[1]]))){
        # somehow we can not loop over the transitions so we use
        # seq_along(length(transitions))
        transition = transitions[[1]][[i]]
        outputs = xmlChildren(transition)$listOfOutputs
        inputs = xmlChildren(transition)$listOfInputs
        signs = as.vector(unlist(sapply(xmlChildren(inputs), function(x) xmlGetAttr(x,"sign"))))
        LHS = as.vector(unlist(sapply(xmlChildren(inputs), function(x) xmlGetAttr(x,"qualitativeSpecies"))))
        RHS = as.vector(unlist(sapply(xmlChildren(outputs), function(x) xmlGetAttr(x,"qualitativeSpecies"))))
        functions = xmlChildren(transition)$listOfFunctionTerms

        # from function, figure out if we have a OR or AND gate
        logics = xmlChildren(xmlChildren(xmlChildren(xmlChildren(functions)$functionTerm)$math)$apply)
        logics = names(logics)

        # all links other than AND should be here
        if (("and" %in% logics)==FALSE){
            for (j in seq_along(LHS)){
                if (signs[[j]] == "positive"){
                    sif = rbind(c(LHS[[j]],1,RHS), sif)
                } else if (signs[[j]] == "negative") {
                    sif = rbind(c(LHS[[j]],-1,RHS), sif)
                } else {
                    stop("signs must be negative or positive")
                }
            }
        } 

        # the AND case only
        if ("and" %in% logics){
            andName = paste("and", andCounter, sep="")
            for (j in seq_along(LHS)){
                if (signs[[j]] == "positive"){
                    sif = rbind(c(LHS[[j]],1,andName), sif)
                } else if (signs[[j]] == "negative") {
                    sif = rbind(c(LHS[[j]],-1,andName), sif)
                } else {
                    stop("signs must be negative or positive")
                }
            }
            sif = rbind(c(andName, 1, RHS), sif)
            andCounter = andCounter + 1
        }
    }

    # remove last row (0,0,0) set at the beginning to define the data frame.
    sif = sif[-dim(sif)[[1]],]

    fh = tempfile()
    write.table(sif,file=fh,
            row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")

    sif = readSIF(fh)

    # hack/trick to replace self loop with a dummy node.
    #sif = .add_dummies(sif)

    return(sif)

}





# this function reads a model, find self loops.
# If any self loops are found, there are replaced by A->dummy->A
# Does handle simple direct self loops e.g. A->A but not A-|A
.add_dummies <- function(model){

    # identify self-loops 
    indices = c()
    for (ireac in seq_along(model$reacID)){
        reac = model$reacID[ireac]
        lhs = strsplit(reac, "=")[[1]][1]
        rhs = strsplit(reac, "=")[[1]][2]
        if (lhs == rhs){
            indices = c(indices, ireac)
        }
    }

    # don't do anything if no self loops are found
    if (length(indices)==0){
        return(model)
    }

    # remove self loops
    warning("Self inhibited loop not handled yet in this version")
    newmodel = model
    newmodel$reacID = newmodel$reacID[-indices]
    newmodel$interMat <- newmodel$interMat[,-indices]
    newmodel$notMat <- newmodel$notMat[,-indices]

    # add them back with the dummies
    counter = 1
    for (i in seq_along(indices)){
        index = indices[i]
        dummy = paste("dummy" , counter, sep="")
        reac = model$reacID[index]
        species = strsplit(reac, "=")[[1]][1]
        reac1 = paste(species, "=", dummy, sep="")
        reac2 = paste(dummy, "=", species, sep="")
        counter = counter + 1

        # add the 2 new reactions only in reacID for now
        newmodel$reacID = c(newmodel$reacID, c(reac1, reac2))

        # add the dummy species
        newmodel$namesSpecies = c(newmodel$namesSpecies, dummy)

        ########## INTERMAT ############
        # now add first the new species row in interMat
        rowvec = rep(0, length(newmodel$reacID)-2 )
        newmodel$interMat = rbind(newmodel$interMat, rowvec)
        names = rownames(newmodel$interMat)
        names[length(names)] = dummy
        rownames(newmodel$interMat) = names

        # add species=dummy column in interMat
        colvec = rep(0, length(newmodel$namesSpecies)) # !! the dummy node is there
        colvec[newmodel$namesSpecies == species] = -1
        colvec[newmodel$namesSpecies == dummy] = 1
        colvec = colvec[1:length(colvec)]
        newmodel$interMat = cbind(newmodel$interMat, colvec)
        names = colnames(newmodel$interMat)
        names[length(names)] = reac1
        colnames(newmodel$interMat) = names

        # add dummy=species column in interMat
        colvec = rep(0, length(newmodel$namesSpecies))
        colvec[newmodel$namesSpecies == species] = 1
        colvec[newmodel$namesSpecies == dummy] = -1
        newmodel$interMat = cbind(newmodel$interMat, colvec)
        names = colnames(newmodel$interMat)
        names[length(names)] = reac2
        colnames(newmodel$interMat) = names


        ########## NOTMAT ############
        # now add the new species row in interMat
        rowvec = rep(0, length(newmodel$reacID)-2)
        newmodel$notMat = rbind(newmodel$notMat, rowvec)
        names = rownames(newmodel$notMat)
        names[length(names)] = dummy
        rownames(newmodel$notMat) = names

        # add species=dummy column in NotMat only zeros
        colvec = rep(0, length(newmodel$namesSpecies))
        newmodel$notMat = cbind(newmodel$notMat, colvec)
        names = colnames(newmodel$notMat)
        names[length(names)] = reac1
        colnames(newmodel$notMat) = names

        ## add dummy=species column in notMat only zeros
        colvec = rep(0, length(newmodel$namesSpecies))
        newmodel$notMat = cbind(newmodel$notMat, colvec)
        names = colnames(newmodel$notMat)
        names[length(names)] = reac2
        colnames(newmodel$notMat) = names


    }
    return(newmodel)
}


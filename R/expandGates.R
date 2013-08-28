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
# $Id: expandGates.R 3859 2013-07-30 13:28:40Z cokelaer $
expandGates<-function(model, ignoreList=NA,maxInputsPerGate=2){

    Model = model
    # check that Model is a Model list
    if(!is.list(Model)) stop("This function expects as input a Model as output by readSIF")
    if(length(Model) == 4) {
        if(all(names(Model) != c("reacID", "namesSpecies","interMat","notMat"))){
            stop("This function expects as input a Model as output by readSIF")
            }
        }
    if(length(Model) == 5) {
        if(all(names(Model) != c("reacID", "namesSpecies","interMat","notMat","speciesCompressed"))) {
            stop("This function expects as input a Model as output by readSIF")
        }
    }

    SplitANDs <- list(initialReac=c("split1","split2"))
    splitR <- 1

    # split all the ANDs
    # remove any ANDs 2/3 and save >3 to add later
    # +3 won't get added here again but are part of prior knowledge if contained within PKN
    andToAdd = c()
    remove.and = c()
    reacs2Ignore = c()
    initialReacN <- length(Model$reacID)

  #TODO: move into readSIF ?
  if (initialReacN == 1){
        Model$interMat <- as.matrix(Model$interMat)
   }

    # which reactions have ignoreList as output?
    if(!is.na(ignoreList[1])) {
        for(s in 1:initialReacN) {
            if(any(Model$interMat[ignoreList,s] == 1)) {
                reacs2Ignore = c(reacs2Ignore, s)
            }
        }
    }

    for(r in 1:initialReacN) {
        inNodes <- which(Model$interMat[,r] == -1)
        if(length(inNodes) > 1) {
            if(length(inNodes) > 3) {
                andToAdd = c(andToAdd, r)
            }
            remove.and = c(remove.and, r)

            if(!any(reacs2Ignore == r)) {
                outNode <- which(Model$interMat[,r] == 1)
                newReacs <- matrix(data=0,nrow=dim(Model$interMat)[1],ncol=length(inNodes))
                newReacsNots <- matrix(data=0,nrow=dim(Model$interMat)[1],ncol=length(inNodes))
                newReacs[outNode,] <- 1
                newReacsIDs <- rep("a",length(inNodes))
                for(i in 1:length(inNodes)) {
                    newReacs[inNodes[i],i] <- -1
                    newReacsNots[inNodes[i],i] <- Model$notMat[inNodes[i],r]
                    newReacsIDs[i] <- paste(Model$namesSpecies[inNodes[i]],"=",Model$namesSpecies[outNode],sep="")
                    if(Model$notMat[inNodes[i],r] == 1) newReacsIDs[i] <- paste("!",newReacsIDs[i],sep="")
                }
                colnames(newReacs) <- newReacsIDs
                colnames(newReacsNots) <- newReacsIDs
                SplitANDs[[splitR]] <- newReacsIDs
                names(SplitANDs)[splitR] <- Model$reacID[r]
                splitR <- splitR+1
                Model$notMat <- cbind(Model$notMat,newReacsNots)
                Model$interMat <- cbind(Model$interMat,newReacs)
                Model$reacID <- c(Model$reacID,newReacsIDs)
            }
        }
    }

    # ***** ADDED *****
    # remove 'AND' gates that will be made anyway
    # save anything > 3 to add at end

    if(length(andToAdd)) {
        toAdd = list()
        toAdd$notMat <- Model$notMat[,andToAdd]
        toAdd$interMat <- Model$interMat[,andToAdd]
        toAdd$reacID <- Model$reacID[andToAdd]
    } else {
        toAdd <- NA
    }

    if(length(remove.and)) {
        Model$notMat <- Model$notMat[,-remove.and]
        Model$interMat <- Model$interMat[,-remove.and]
        Model$reacID <- Model$reacID[-remove.and]
    }

    # ***** ADDED *****

    # the list SplitANDs now contains a named element for each AND reac that has been split,
    # and each element contains a vector with the names of the reactions that result from the split
    # create combinations of ORs
    # the newANDs list will contain an element for each new '&' gate, named by the name of this new and
    # reac, and containing a vector of the names of the reactions from which it was created

    newANDs <- list(finalReac=c("or1","or2"))
    ANDsadded <- 1
    total.list = 1:length(Model$namesSpecies)


    # functions to get lhs and rhs of reactions
    getlhs <- function(x) {
        spec1 = strsplit(x, "=")[[1]][1]
    }
    getrhs <- function(x) {
        spec2 = strsplit(x, "=")[[1]][2]
    }

    # scan all species and build and gates if required
    for(sp in total.list) {
        inReacsIndex <- which(Model$interMat[sp,] == 1)
        if(length(inReacsIndex) > 1) {
            inReacs <- Model$interMat[,inReacsIndex]
            findInput <- function(x) {
                inNode<-which(x == -1)
            }
            inSp <- apply(inReacs,2,findInput)

            # let
            # find the input species first and store in a vector
            inSpecies = apply(as.matrix(colnames(inReacs)), 1, getlhs)
            outname = Model$namesSpecies[sp]

            # just for sanity check, all outputs must be the same
            outnames = apply(as.matrix(colnames(inReacs)), 1, getrhs)
            if (length(unique(outnames))!=1 | outname!=outnames[1]){
                stop("error in expandGates. should not happen here.please report")
            }

            # an alias
            myrownames = rownames(Model$interMat)

            # first the 2 inputs cases

            combinations = combn(seq(1,length(inSpecies)), 2)

            for (this in seq(1, dim(combinations)[2])){
                i = combinations[1,this]
                j = combinations[2,this]
                # names[i] and names[j] contains possibly the !
                # sign,let us get the real species names
                realname1 = ifelse(substr(inSpecies[i], 1,1)=="!",substr(inSpecies[i],2,10000), inSpecies[i])
                realname2 = ifelse(substr(inSpecies[j], 1,1) =="!",substr(inSpecies[j],2,10000), inSpecies[j])

                realnames = c(realname1,realname2)
                if (any(combn(realnames,2)[1,] == combn(realnames,2)[2,])){  # exclude reaction if any name are indentical
                    next()
                }

                # create the new reaction Id to be used as a column name
                newcolname = paste(paste(inSpecies[i], inSpecies[j],sep="+"),outname, sep="=")
                if (newcolname %in% colnames(Model$interMat)){
                    next() # skip if exist already
                }
                Model$reacID <- c(Model$reacID,newcolname)

                # fill the interMat (-1 if in lhs, 1 if in rhs)
                values = as.matrix(rep(0, length(Model$namesSpecies)))
                colnames(values)<-newcolname
                values[which(myrownames == realname1)]<- -1
                values[which(myrownames == realname2)]<- -1
                values[which(myrownames == outname)]<- 1
                Model$interMat= cbind(Model$interMat, values)

                # now, the notMat, 0 by default
                values = as.matrix(rep(0, length(Model$namesSpecies)))
                colnames(values)<-newcolname
                if (substr(inSpecies[i],1,1) == "!"){ #look first specy
                    values[which(myrownames == realname1)]<- 1
                }
                if (substr(inSpecies[j],1,1) == "!"){# and second one
                    values[which(myrownames == realname2)]<- 1
                }
                Model$notMat= cbind(Model$notMat, values)

                # finally, fill the newAnd list
                newreac1 = paste(inSpecies[i], outname, sep="=")
                newreac2 = paste(inSpecies[j], outname, sep="=")
                newANDs[[length(newANDs)+1]] <- c(newreac1, newreac2)
                names(newANDs)[[length(newANDs)]] <- newcolname
            }

            # Same code as above but to create the 3 inputs combinations
            if (length(inSpecies)>=3 & maxInputsPerGate>=3){
                combinations = combn(seq(1,length(inSpecies)), 3)
                indices = seq(1,dim(combinations)[2])
            }
            else{
                indices = seq(length=0)
            }

            for (this in indices){
                i = combinations[1,this]
                j = combinations[2,this]
                k = combinations[3,this]
                realname1 = ifelse(substr(inSpecies[i], 1,1)=="!",substr(inSpecies[i],2,10000), inSpecies[i])
                realname2 = ifelse(substr(inSpecies[j], 1,1) =="!",substr(inSpecies[j],2,10000), inSpecies[j])
                realname3 = ifelse(substr(inSpecies[k], 1,1) =="!",substr(inSpecies[k],2,10000), inSpecies[k])

                realnames = c(realname1,realname2, realname3)
                if (any(combn(realnames,2)[1,] == combn(realnames,2)[2,])){  # exclude reaction if any name are indentical
                    next()
                }
                newcolname <- paste(paste(inSpecies[i], inSpecies[j],inSpecies[k],
                               sep="+"), outname, sep="=")
                if (newcolname %in% colnames(Model$interMat)){
                    next() # skip if exist already
                }
                Model$reacID <- c(Model$reacID,newcolname)

                # intermat first
                values = as.matrix(rep(0, length(Model$namesSpecies)))
                colnames(values)<-newcolname
                for (name in inSpecies){
                    realname = ifelse(substr(name, 1,1)=="!",substr(name,2,10000), name)
                    values[which(myrownames == realname)]<- -1
                }
                values[which(myrownames == outname)]<- 1
                Model$interMat= cbind(Model$interMat, values)

                # now, the notMat
                values = as.matrix(rep(0, length(Model$namesSpecies)))
                colnames(values)<-newcolname
                for (name in inSpecies){
                    if (substr(name,1,1) == "!"){
                        realname = ifelse(substr(name, 1,1) =="!", substr(name,2,10000), name)
                        values[which(myrownames == realname)]<- 1
                    }
                }
                Model$notMat= cbind(Model$notMat, values)

                # finally the newAnd
                newreac1 = paste(inSpecies[i], outname, sep="=")
                newreac2 = paste(inSpecies[j], outname, sep="=")
                newreac3 = paste(inSpecies[k], outname, sep="=")
                newANDs[[length(newANDs)+1]] <- c(newreac1, newreac2, newreac3)
                names(newANDs)[[length(newANDs)]] <- newcolname
            }

            # Same code as above but to create the 4 inputs combinations
            if (length(inSpecies)>=4 & maxInputsPerGate>=4){
                combinations = combn(seq(1,length(inSpecies)), 4)
                indices = seq(1,dim(combinations)[2])
            }
            else{
                indices = seq(length=0)
            }

            for (this in indices){
                i = combinations[1,this]
                j = combinations[2,this]
                k = combinations[3,this]
                l = combinations[4,this]
                realname1 = ifelse(substr(inSpecies[i], 1,1)=="!",substr(inSpecies[i],2,10000), inSpecies[i])
                realname2 = ifelse(substr(inSpecies[j], 1,1) =="!",substr(inSpecies[j],2,10000), inSpecies[j])
                realname3 = ifelse(substr(inSpecies[k], 1,1) =="!",substr(inSpecies[k],2,10000), inSpecies[k])
                realname4 = ifelse(substr(inSpecies[l], 1,1) =="!",substr(inSpecies[l],2,10000), inSpecies[l])
                realnames = c(realname1,realname2, realname3, realname4)
                if (any(combn(realnames,2)[1,] == combn(realnames,2)[2,])){  # exclude reaction if any name are indentical
                    next()
                }
                newcolname <- paste(paste(inSpecies[i], inSpecies[j],inSpecies[k],inSpecies[l],
                               sep="+"), outname, sep="=")
                if (newcolname %in% colnames(Model$interMat)){
                    next() # skip if exist already
                }
                Model$reacID <- c(Model$reacID,newcolname)

                # intermat first
                values = as.matrix(rep(0, length(Model$namesSpecies)))
                colnames(values)<-newcolname
                for (name in inSpecies){
                    realname = ifelse(substr(name, 1,1)=="!",substr(name,2,10000), name)
                    values[which(myrownames == realname)]<- -1
                }
                values[which(myrownames == outname)]<- 1
                Model$interMat= cbind(Model$interMat, values)

                # now, the notMat
                values = as.matrix(rep(0, length(Model$namesSpecies)))
                colnames(values)<-newcolname
                for (name in inSpecies){
                    if (substr(name,1,1) == "!"){
                        realname = ifelse(substr(name, 1,1) =="!", substr(name,2,10000), name)
                        values[which(myrownames == realname)]<- 1
                    }
                }
                Model$notMat= cbind(Model$notMat, values)

                # finally the newAnd
                newreac1 = paste(inSpecies[i], outname, sep="=")
                newreac2 = paste(inSpecies[j], outname, sep="=")
                newreac3 = paste(inSpecies[k], outname, sep="=")
                newreac4 = paste(inSpecies[l], outname, sep="=")
                newANDs[[length(newANDs)+1]] <- c(newreac1, newreac2, newreac3, newreac4)
                names(newANDs)[[length(newANDs)]] <- newcolname

            } # end if length(inSp) == 2


 if (maxInputsPerGate >= 5) {
        for (mip in 5:maxInputsPerGate) {
          if (length(inSpecies) >= mip & maxInputsPerGate >= 
              mip) {
            combinations = combn(seq(1, length(inSpecies)), 
              mip)
            indices = seq(1, dim(combinations)[2])
          }
          else {
            indices = seq(length = 0)
          }
          for (this in indices) {
            combs <- c()
            realnames <- c()
            for (i in 1:mip) {
              combs[i] = combinations[i, this]
              realnames[i] = ifelse(substr(inSpecies[combs[i]], 1, 1) == 
                         "!", substr(inSpecies[i], 2, 10000), inSpecies[combs[i]])
            }
            if (any(combn(realnames, 2)[1, ] == combn(realnames, 
                                         2)[2, ])) {
              (next)()
            }
            newcolname <- paste(paste(inSpecies[combs], collapse = "+"), outname, 
                                sep = "=")
            if (newcolname %in% colnames(Model$interMat)) {
              (next)()
            }
            Model$reacID <- c(Model$reacID, newcolname)
            values = as.matrix(rep(0, length(Model$namesSpecies)))
            colnames(values) <- newcolname
            for (name in inSpecies) {
              realname = ifelse(substr(name, 1, 1) == "!", 
                substr(name, 2, 10000), name)
              values[which(myrownames == realname)] <- -1
            }
            values[which(myrownames == outname)] <- 1
            Model$interMat = cbind(Model$interMat, values)
            values = as.matrix(rep(0, length(Model$namesSpecies)))
            colnames(values) <- newcolname
            for (name in inSpecies) {
              if (substr(name, 1, 1) == "!") {
                realname = ifelse(substr(name, 1, 1) == "!", 
                  substr(name, 2, 10000), name)
                values[which(myrownames == realname)] <- 1
              }
            }
            Model$notMat = cbind(Model$notMat, values)
            newreac1 = paste(inSpecies[i], outname, sep = "=")
            newreac2 = paste(inSpecies[j], outname, sep = "=")
            newreac3 = paste(inSpecies[k], outname, sep = "=")
            newreac4 = paste(inSpecies[l], outname, sep = "=")
            newreacs <- c()
            for (i in 1:mip) {
              newreacs[i] <- paste(inSpecies[combs[i]], outname, sep = "=")
            }
            newANDs[[length(newANDs) + 1]] <- newreacs # c(newreac1, newreac2, newreac3, newreac4)
            names(newANDs)[[length(newANDs)]] <- newcolname
          }
        }
      }

        }
    }
#    dup.index2 = c(1:length(Model$reacID))
#    dup.index2 = dup.index2[duplicated(Model$reacID)]
#    Model$reacID = Model$reacID[-dup.index2]
#    Model$interMat = Model$interMat[,-dup.index2]
#    Model$notMat = Model$notMat[,-dup.index2]

    if(!is.na(toAdd)) {
        Model$notMat = cbind(Model$notMat, toAdd$notMat)
        Model$interMat = cbind(Model$interMat, toAdd$interMat)
        Model$reacID = c(Model$reacID, toAdd$reacID)
    }

    modelExp <- Model
    modelExp$SplitANDs <- SplitANDs
    modelExp$newANDs <- newANDs
    return(modelExp)
}


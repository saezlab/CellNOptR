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
# $Id: compressModel.R 590 2012-02-22 17:16:50Z cokelaer $
compressModel<-function(Model, indexes){

    #this will hold the indices of the species that we have to look at (i.e. not signal nor cue)
    species<-seq(1:length(Model$namesSpecies))
    #species<-species[-c(indexes$signals,indexes$stimulated,indexes$inhibited)]
    # TODO: is.na proposed by Aidan fix. 
    if(is.na(indexes$inhibited[1])) {
        species <- species[-c(indexes$signals,indexes$stimulated)]
    } else {
        species <- species[-c(indexes$signals,indexes$stimulated,indexes$inhibited)]
    }

    speciesCompressed<-0

    for(sp in species){
        # Find all the reacs that come out of species sp
        outReac <- which(Model$interMat[sp,] == -1)

        # Find all the reacs that come into species sp
        inReac <- which(Model$interMat[sp,] == 1)

        # If the species doesn't have any output, we can't compress it, so we do sthg only if
        #the species has one output minimum
        if(length(outReac) != 0){

            #If the species has only one input        
            if(length(inReac) == 1 && length(which(Model$interMat[,inReac] == -1)) == 1) { 
                #TODO aidan fix why was the second condition added?

				# find the input element
                input <- which(Model$interMat[,inReac] == -1)

                # find the output elements
                if(length(outReac) > 1) outputs <- apply(Model$interMat[,outReac],2,function(x) which(x == 1))
                if(length(outReac) == 1) outputs <- which(Model$interMat[,outReac] == 1)
            
                # check that the input is not an & gate, and that none of the downstream
                # elements is the output element
                if(length(input) == 1 && !any(outputs == input[1]) 
                    && !any(abs(Model$interMat[,inReac] + Model$interMat[,outReac[1]]) > 1)) { 
                    # TODO: fix properly: the case where the input is also part of an 'AND' gate in the output. outcome: don't compress
                    # need to create new reacs, that go from input to outputs, and to remove the old reacs

                    # Create the first one, that goes from input to the first output
                    reac2Remove <- c(inReac, outReac)
                    newReacs = matrix(NA, nrow=dim(Model$interMat)[1], ncol=length(outReac))
                    rownames(newReacs) = Model$namesSpecies
                    newNots = matrix(NA, nrow=dim(Model$notMat)[1], ncol=length(outReac))
                    rownames(newNots) = Model$namesSpecies
                    newReacIDs = rep(NA, length(outReac))
                    
                    for(r in 1:length(outReac)) {
                        newReacs[,r] <- Model$interMat[,inReac] + Model$interMat[,outReac[r]]
                        # TODO: check this?
                        newNots[,r] <- Model$notMat[,outReac[r]]
                        if(Model$notMat[sp,outReac[r]] == 1 && Model$notMat[which(Model$interMat[,inReac]==-1),inReac]==1) 
                        newNots[input,r] <- 0
                        if(Model$notMat[sp,outReac[r]] == 0 && Model$notMat[which(Model$interMat[,inReac]==-1),inReac]==0) 
                        newNots[input,r] <- 0
                        if(Model$notMat[sp,outReac[r]] == 1 && Model$notMat[which(Model$interMat[,inReac]==-1),inReac]==0) 
                        newNots[input,r] <- 1
                        if(Model$notMat[sp,outReac[r]] == 0 && Model$notMat[which(Model$interMat[,inReac]==-1),inReac]==1) 
                        newNots[input,r] <- 1
                        if(length(grep("+",names(outReac[r])))) {
                            andInput = rownames(newReacs)[which(newReacs[,r] == -1)]
                            andInd = which(newReacs[,r] == -1)
                            if(length(intersect(which(newNots[,r]==1),which(newReacs[,r]==-1)))) {
                                andNeg = intersect(which(newNots[,r]==1),which(newReacs[,r]==-1))
                                for(p in 1:length(andNeg)) {
                                    andInput[which(andInd==andNeg[p])] = paste("!", andInput[which(andInd==andNeg[p])], sep="")
                                }
                            }
                            
                            LHS = paste(andInput,collapse="+", sep="")
                            newReacIDs[r] = paste(LHS, "=", rownames(newReacs)[which(newReacs[,r] == 1)], sep="")
                            
                        } else {
                            newReacIDs[r] <- paste(rownames(newReacs)[which(newReacs[,r] == -1)], "=", rownames(newReacs)[which(newReacs[,r] == 1)],sep="")
                            if(any(newNots[,r] == 1)) newReacIDs[r] <- paste("!", newReacIDs[4], sep="")
                        }    

                    }

                    colnames(newReacs) <- newReacIDs
                    colnames(newNots) <- newReacIDs
                    
                    # now remove the reacs to remove and append the ones that have just been created
                    Model$interMat <- Model$interMat[,-reac2Remove]
                    Model$notMat <- Model$notMat[,-reac2Remove]
                    Model$reacID <- Model$reacID[-reac2Remove]
                    Model$interMat <- cbind(Model$interMat, newReacs)
                    Model$notMat <- cbind(Model$notMat,newNots)
                    Model$reacID <- c(Model$reacID,newReacIDs)
                
                    if(length(outReac) == 1) {
                        colnames(Model$interMat)[dim(Model$interMat)[2]] <- newReacIDs
                        colnames(Model$notMat)[dim(Model$notMat)[2]] <- newReacIDs
                    }

                    # store which species has been compressed. This will be used to remove rows of the matrices
                    # after all is compressed (I don't do it know because then the indexes of the species that we
                    # have decided to look at in the beginning won't be valid anymore)
                
                    speciesCompressed <- c(speciesCompressed,sp)
                } #end of if(length(input) == 1 && !any(outputs == input[1]) && !any....
            } else {

                # If the species has more than one input, but only one output
                if(length(outReac) == 1 && length(which(Model$interMat[,outReac] == -1)) == 1) {
                    
                    # find the output element
                    output <- which(Model$interMat[,outReac] == 1)
                    # find the input element(s)
                    inputs = list()
                    
                    if(length(inReac) > 1) {
                        for(a in inReac) {
                            inputs = c(inputs, list(which(Model$interMat[,a] == -1)))
                        }
                    } 
                    if(length(inReac) == 1) inputs <- list(which(Model$interMat[,inReac] == -1))

                    # check that none of the downstream elements is the output element
                    inputsUnlist = unlist(inputs)
                    if(!any(inputsUnlist == output[1])) {
                        # create the first reac, that goes from the first input to the output
                        reac2Remove <- c(inReac, outReac)
                        newReacs = matrix(NA, nrow=dim(Model$interMat)[1], ncol=length(inReac))
                        rownames(newReacs) = Model$namesSpecies
                        newNots = matrix(NA, nrow=dim(Model$notMat)[1], ncol=length(inReac))
                        rownames(newNots) = Model$namesSpecies
                        newReacIDs = rep(NA, length(inReac))
                        if (length(inReac)>0){
                        for(r in 1:length(inReac)) {
                            newReacs[,r] <- Model$interMat[,inReac[r]] + Model$interMat[,outReac]
                            newNots[,r] <- Model$notMat[,inReac[r]]
                            # if out is negative & in is single:
                            if(length(inputs[[r]])==1 && Model$notMat[sp,outReac]==1) {
                                newNots[inputs[[r]],r] <- ifelse((newNots[inputs[[r]],r] == 1), 0, 1)
                            }
                            # if out is negative & in is 'AND' gate
                            if(length(inputs[[r]]) > 1 && Model$notMat[sp,outReac]==1) {
                                for(b in 1:length(inputs[[r]])) {
                                    newNots[inputs[[r]][b]] <- ifelse((newNots[inputs[[r]][b]] == 1), 0, 1)
                                }
                            }
                                    
                            if(length(grep("+",names(inReac[r])))) {
                                andInput = rownames(newReacs)[which(newReacs[,r] == -1)]
                                andInd = which(newReacs[,r] == -1)
                                if(length(intersect(which(newNots[,r]==1),which(newReacs[,r]==-1)))) {
                                    andNeg = intersect(which(newNots[,r]==1),which(newReacs[,r]==-1))
                                    for(p in 1:length(andNeg)) {
                                        andInput[which(andInd==andNeg[p])] = paste("!", andInput[which(andInd==andNeg[p])], sep="")
                                    }
                                }
                            
                            LHS = paste(andInput,collapse="+",sep="")
                            newReacIDs[r] = paste(LHS, "=", rownames(newReacs)[which(newReacs[,r] == 1)], sep="")
                            
                            } else {
                                newReacIDs[r] <- paste(rownames(newReacs)[which(newReacs[,r] == -1)], "=", rownames(newReacs)[which(newReacs[,r] == 1)],sep="")
                                if(any(newNots[,r] == 1)) newReacIDs[r] <- paste("!", newReacIDs[4], sep="")
                            }    

                        }}
                        
                        colnames(newReacs) <- newReacIDs
                        colnames(newNots) <- newReacIDs
                        
                        Model$interMat <- Model$interMat[,-reac2Remove]
                        Model$notMat <- Model$notMat[,-reac2Remove]
                        Model$reacID <- Model$reacID[-reac2Remove]
                        Model$interMat <- cbind(Model$interMat, newReacs)
                        Model$notMat <- cbind(Model$notMat, newNots)
                        Model$reacID <- c(Model$reacID, newReacIDs)
                        if(length(inReac) == 1) {
                            colnames(Model$interMat)[dim(Model$interMat)[2]] <- newReacIDs
                            colnames(Model$notMat)[dim(Model$notMat)[2]] <- newReacIDs
                        }    
                
                        speciesCompressed <- c(speciesCompressed,sp)    
                    }
                }
            }

        # If the species has more than one input and more than one output then
        # nothing is done.
        }
    } # end of the for loop over species



    ModelCompressed<-Model
    ModelCompressed$speciesCompressed<-NA

    # If species have been compressed
    if(length(speciesCompressed) > 1){

        # Remove the dummy 0 in speciesCompressed
        speciesCompressed<-speciesCompressed[2:length(speciesCompressed)]
        ModelCompressed$speciesCompressed<-Model$namesSpecies[speciesCompressed]
        ModelCompressed$namesSpecies<-Model$namesSpecies[-speciesCompressed]
        ModelCompressed$notMat<-Model$notMat[-speciesCompressed,]
        ModelCompressed$interMat<-Model$interMat[-speciesCompressed,]
    }

    # Check that no duplicates reactions have been created
    uniqueR<-unique(ModelCompressed$reacID)

    if(length(uniqueR) != length(ModelCompressed$reacID)){
        ModelCompressed$reacID<-ModelCompressed$reacID[match(uniqueR,ModelCompressed$reacID)]
        ModelCompressed$interMat<-ModelCompressed$interMat[,match(uniqueR,colnames(ModelCompressed$interMat))]
        ModelCompressed$notMat<-ModelCompressed$notMat[,match(uniqueR,colnames(ModelCompressed$notMat))]
        }

    return(ModelCompressed)    

}


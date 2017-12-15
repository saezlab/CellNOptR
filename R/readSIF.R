# -*- python -*-
#
#  This file is part of the CNO software
#
#  Copyright (c) 2011-2012 - EMBL - European Bioinformatics Institute
#
#  File author(s): Thomas Cokelaer <cokelaer@ebi.ac.uk>
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  CNO website: http://www.cellnopt.org
#
##############################################################################
# $Id$
readSIF<-function(sifFile){
	
    # Read the sif file expected 3 columns all over the file
    status = tryCatch({
            sif<-read.table(sifFile)
            sif<-as.matrix(sif)
            status = TRUE},
        error = function(e) { status = FALSE })

    # if the sif file has more than 3 columns e.g.:
    # A 1 B C
    # B 1 D
    # then another parser is required:
    if (status == FALSE){
        data = readLines(sifFile)
        # replace all tabs by spaces
        data = gsub("\t", " ",data)
        # replace all long whitespaces by a single space character
        data = gsub("  *", " ",data)
        sif = c()
        # scan all lines and extract the reactions
        for (line in data){
            items = unlist(strsplit(line, " "))
            LHS = items[1]
            edge = items[2]
            RHSs = items[-c(1,2)] # let us look at the RHS elements 
            for (RHS in RHSs){
                sif = c(sif, LHS)
                sif = c(sif, edge)
                sif = c(sif, RHS)
            }
        }
        sif = matrix(sif, ncol=3, byrow=TRUE)
    }
	
	sif <- unique(sif)   # we only select the unique interactions

    #Create the vector of names of species
    namesSpecies<-unique(c(as.character(sif[,1]), as.character(sif[,3])))

    #Create the vector of names or reactions

    createReacID<-function(x){
        neg<-ifelse(test=(x[2] == "-1"), TRUE, FALSE)
        ID<-ifelse(
            neg,
            paste("!",x[1],"=",x[3],sep=""),
            paste(x[1],"=",x[3],sep=""))
        return(ID)
        }

    reacID<-apply(sif,1,createReacID)

    # reaction with and gates can not have "and -1 specy" reaction: edge must be
    # 1
    if (any(substring(reacID, 1,4)=="!and")==TRUE){
        stop("Found an and gate with not edge. This is not possible check your SIF file. Must be 'and 1 specy', not 'and -1 specy' !!")
    }
    if (any(substring(reacID, 1,4)=="!AND")==TRUE){
        stop("Found an and gate with not edge. This is not possible check your SIF file. Must be 'and 1 specy', not 'and -1 specy' !!")
    }
#Create empty interMat and notMat
    interMat<-matrix(data=0,nrow=length(namesSpecies), ncol=dim(sif)[1])
    rownames(interMat)<-namesSpecies
    colnames(interMat)<-reacID
    notMat<-matrix(data=0,nrow=length(namesSpecies), ncol=dim(sif)[1])
    rownames(notMat)<-namesSpecies
    colnames(notMat)<-reacID


#Fill in interMat and notMat
    for(r in 1:length(reacID)){

        source<-match(sif[r,1],namesSpecies)
        target<-match(sif[r,3],namesSpecies)
        neg<-ifelse(test=(sif[r,2] == "-1"), TRUE, FALSE)
        interMat[source,r]<-(-1)
        interMat[target,r]<-1
        if(neg) notMat[source,r]<-1

        }

#Take care of the "and" nodes
#browser()
#1.detect them
    AndNodes<-grep(
        pattern="(^[a,A][n,N][d,D]\\d+$)",namesSpecies,perl=TRUE,ignore.case=FALSE)
    AndNodesV<-grep(
        pattern="(^[a,A][n,N][d,D]\\d+$)",namesSpecies,perl=TRUE,ignore.case=FALSE,value=TRUE)


    if(length(AndNodes) != 0){

#2.go through them
        for(i in 1:length(AndNodes)){

#3.store the name of the node, and find the associated reactions
            node<-AndNodesV[i]
            AndsReacs<-c(
                grep(pattern=paste("^",node,"=",sep=""),reacID,ignore.case=FALSE),
                grep(pattern=paste(node,"$",sep=""),perl=TRUE,reacID,ignore.case=FALSE))

#4.create the new column of interMat and notMat
            newReac<-rowSums(interMat[,AndsReacs])
            newReacNot<-rowSums(notMat[,AndsReacs])

#5.create the new reacID

            #First case: I have an "and" with two inputs only

            if(length(AndsReacs) == 3){

                output<-grep(pattern=paste("^",node,"=",sep=""),reacID,
                    ignore.case=FALSE,value=TRUE)
                inputs<-grep(pattern=paste(node,"$",sep=""),perl=TRUE,
                    reacID,ignore.case=FALSE,value=TRUE)
                output<-sub(pattern=node,x=output,replacement="",ignore.case=FALSE)
                inputs<-sub(pattern=paste("=",node,"$",sep=""),x=inputs,
                    replacement="+",ignore.case=FALSE)
                inputs<-paste(inputs,collapse='')
                newReacID<-paste(inputs,output,sep="")
                newReacID<-sub(pattern=paste("\\W{1}",output,sep=""),
                    x=newReacID,replacement=output,ignore.case=FALSE,perl=TRUE)

                }else{

            #Second case: I have an 'and' with more than 2 inputs
                    output<-grep(pattern=paste("^",node,"=",sep=""),reacID,
                        ignore.case=FALSE,value=TRUE)

                    if(length(output) == 1){

                        output<-sub(pattern=node,x=output,replacement="",ignore.case=FALSE)
                        inputs<-grep(pattern=paste(node,"$",sep=""),perl=TRUE,
                            reacID,ignore.case=FALSE,value=TRUE)
                        inputs<-sub(pattern=paste("=",node,"$",sep=""),x=inputs,
                            replacement="+",ignore.case=FALSE)
                        inputs<-paste(inputs,collapse='')
                        newReacID<-paste(inputs,output,sep="")
                        newReacID<-sub(pattern=paste("\\W{1}",output,sep=""),
                            x=newReacID,replacement=output,ignore.case=FALSE,perl=TRUE)

                        }else{

                            output<-sub(pattern=paste("^",node,"=",sep=""),x=output,
                                replacement="",ignore.case=FALSE)
                            inputs<-grep(pattern=paste(node,"$",sep=""),perl=TRUE,
                                reacID,ignore.case=FALSE,value=TRUE)
                            inputs<-sub(pattern=paste("=",node,"$",sep=""),x=inputs,
                                replacement="+",ignore.case=FALSE)
                            inputs<-paste(inputs,collapse='')
                            inputs<-sub(pattern="\\W{1}$",x=inputs,replacement="",
                                ignore.case=FALSE,perl=TRUE)
                            newReacID<-paste(inputs,output,sep="=")

                            }
                    }

#5b.if there are more than one outputs, then the interMat and notMat created above are not correct
            output<-grep(pattern=paste("^",node,"=",sep=""),reacID,
                ignore.case=FALSE,value=FALSE)
            inputs<-grep(pattern=paste(node,"$",sep=""),perl=TRUE,
                reacID,ignore.case=FALSE,value=FALSE)

            if(length(output) != 1){

                newReac<-rowSums(interMat[,c(inputs,output[1])])
                newReacNot<-rowSums(notMat[,c(inputs,output[1])])

                for(n in 2:length(output)){
                    newReac<-cbind(newReac,rowSums(interMat[,c(inputs,output[n])]))
                    newReacNot<-cbind(newReacNot,rowSums(notMat[,c(inputs,output[n])]))
                    }

                }

#6. remove the old reactions and replace them by the new ones
            namesSpecies<-namesSpecies[-which(namesSpecies == AndNodesV[i])]
            reacID<-reacID[-AndsReacs]
            interMat<-interMat[,-AndsReacs,drop=FALSE]
            notMat<-notMat[,-AndsReacs,drop=FALSE]
            interMat<-cbind(interMat,newReac)
            notMat<-cbind(notMat,newReacNot)

            reacID<-c(reacID,newReacID)
            interMat<- as.matrix(interMat[-which(rownames(interMat) == AndNodesV[i]),], drop=F)
            notMat <- as.matrix(notMat[-which(rownames(notMat) == AndNodesV[i]),], drop=F)

            
            if(length(output) == 1){

                colnames(interMat)[dim(interMat)[2]]<-newReacID
                colnames(notMat)[dim(notMat)[2]]<-newReacID

                }else{

                    colnames(interMat)[(dim(interMat)[2]-length(newReacID)+1):dim(interMat)[2]]<-newReacID
                    colnames(notMat)[(dim(notMat)[2]-length(newReacID)+1):dim(notMat)[2]]<-newReacID

                }
            }
        }

    model<-list(
        reacID=reacID,
        namesSpecies=namesSpecies,
        interMat=interMat,
        notMat=notMat)
    return(model)
    }


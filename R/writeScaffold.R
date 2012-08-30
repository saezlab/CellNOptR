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
# $Id: writeScaffold.R 2267 2012-08-30 15:31:54Z cokelaer $
#the global function will still be called writeScaffold and will still need the same arguments
#(ModelComprExpanded,optimResT1,optimResT2,ModelOriginal=ToyModel,CNOlist), but it will be divided into a function writeScaffoldW that does
#the actual writing to files,
#a function that gets the info for the sif file,
#and a function that gets the additional info for the dot file
writeScaffold<-function(
    modelComprExpanded,
    optimResT1,
    optimResT2,
    modelOriginal,
    CNOlist,
    tag=NULL){

#get the stuff that I need for the sif file
    sif<-getSifInfo(modelComprExpanded=modelComprExpanded,
        optimResT1=optimResT1,
        optimResT2=optimResT2,
        modelOriginal=modelOriginal)

#get the stuff that I need for the dot file
    dot<-getDotInfo(
        modelComprExpanded=modelComprExpanded,
        modelOriginal=modelOriginal,
        CNOlist=CNOlist,
        sifFile=sif$sifFile)

#this bit writes the dot, the sif and the sif edges attributes
    writeScaffoldW(
        dN=dot$dN,
        dM=dot$dM,
        modelComprExpanded=modelComprExpanded,
        sifFile=sif$sifFile,
        EApresent=sif$EApresent,
        EAweights=sif$EAweights,
        tag=tag)

    }


###########################################################################
#######these are the functions used above
###########################################################################
#this function writes the sif file, the edge attributes, and the dot file
writeScaffoldW<-function(
    dN,
    dM,
    modelComprExpanded,
    sifFile,
    EApresent,
    EAweights,
    tag=NULL){

    create_filename<-function(x, tag=NULL){
        if (is.null(tag)){
            return(x)
        }
        else{
            return(paste(c(tag, "_", x), collapse=""))
        }
    }

    writeDot(
        dotNodes=dN,
        dotMatrix=dM,
        model=modelComprExpanded,
        filename=create_filename("Scaffold.dot", tag=tag))

    write.table(
        sifFile[,1:3],
        file=create_filename("Scaffold.sif", tag=tag),
        row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")

    write.table(
        EApresent,
        file=create_filename("TimesScaffold.EA", tag=tag),
        row.names=FALSE,col.names="Times",quote=FALSE,sep="\t")

    write.table(
        EAweights,
        file=create_filename("weightsScaffold.EA", tag=tag),
        row.names=FALSE,col.names="Weights",quote=FALSE,sep="\t")

    }

######
#this function computes the stuff that is needed for the dot file

getDotInfo<-function(modelComprExpanded,modelOriginal,CNOlist,sifFile){

    dM<-sifFile
    nodesCompr<-modelComprExpanded$speciesCompressed
    indexes<-indexFinder(CNOlist,modelOriginal)
    nodesSig<-modelOriginal$namesSpecies[indexes$signals]
    nodesInh<-modelOriginal$namesSpecies[indexes$inhibited]
    nodesStim<-modelOriginal$namesSpecies[indexes$stimulated]
    nodesNCNO<-findNONC(modelOriginal,indexes)
    nodesAttrNames<-nodesSig
    nodesAttr<-rep("signal",length(nodesSig))

    if(length(nodesInh) != 0){
        nodesAttrNames<-c(nodesAttrNames,nodesInh)
        nodesAttr<-c(nodesAttr,rep("inhibited",length(nodesInh)))
    }

    if(length(nodesStim) != 0){
        nodesAttrNames<-c(nodesAttrNames,nodesStim)
        nodesAttr<-c(nodesAttr,rep("stimulated",length(nodesStim)))
    }

    if(length(nodesNCNO) != 0){
        nodesAttrNames<-c(nodesAttrNames,nodesNCNO)
        nodesAttr<-c(nodesAttr,rep("ncno",length(nodesNCNO)))
    }

    if(length(nodesCompr) != 0){
        nodesAttrNames<-c(nodesAttrNames,nodesCompr)
        nodesAttr<-c(nodesAttr,rep("compressed",length(nodesCompr)))
    }

    dN<-cbind(nodesAttrNames,nodesAttr)
    return(list(dN=dN,dM=dM))
    }

#####
#this function computes the stuff that is needed for the sif file
getSifInfo<-function(modelComprExpanded,
    optimResT1,
    optimResT2,
    modelOriginal){

    bString1<-optimResT1$bString

    if(is.na(optimResT2[1])){
        bString2<-optimResT1$bString[which(optimResT1$bString == 0)]
        }else{
            bString2<-optimResT2$bString
            }

#BStimes is a string containing 0,1 and 2 depending on whether the interaction is
#absent, present at t1 or present at t2
    BStimes<-bString1
    BStimes[which(BStimes == 0)]<-bString2*2

#create: weightsE is a string that holds the weights of the interactions

    if(!is.null(dim(optimResT1$stringsTol))){
        bW1<-apply(optimResT1$stringsTol,2,mean)
        }else{
            bW1<-bString1
            }

    if(!is.na(optimResT2[1])){

        if(!is.null(dim(optimResT2$stringsTol))){
            bW2<-apply(optimResT2$stringsTol,2,mean)
            }else{
                bW2<-bString2
                }

            weightsE<-bW1
            weightsE[which(optimResT1$bString == 0)]<-weightsE[which(optimResT1$bString == 0)]+bW2

            }else{
                weightsE<-bW1
                }

#These mini functions are used to find the inputs and output of reactions
    findOutput<-function(x){
        sp<-which(x == 1)
        sp<-modelComprExpanded$namesSpecies[sp]
        }

    reacOutput<-apply(modelComprExpanded$interMat,2,findOutput)

    findInput<-function(x){
        sp<-which(x == -1)
        sp<-modelComprExpanded$namesSpecies[sp]
        }

    reacInput<-apply(modelComprExpanded$interMat,2,findInput)

#This mini function is used to create a reaction label as used in a cystoscape edge attribute file
    createReac<-function(x){
        r<-paste(x[1]," (",x[2],") ",x[3],sep="")
        return(r)
        }

#if the class of reacInput is not a list, then there are no AND gates
    if(class(reacInput) != "list"){
        isNeg<-function(x){
            isNegI<-any(x == 1)
            return(isNegI)
        }
        inpSign<-apply(modelComprExpanded$notMat,2,isNeg)
        inpSign<-!inpSign
        inpSign[inpSign]<-1
        inpSign[!inpSign]<--1
        sifFile<-cbind(reacInput,inpSign,reacOutput)
        EApresent<-apply(sifFile,1,createReac)
        EApresent<-cbind(EApresent,BStimes)
        EAweights<-cbind(EApresent,weightsE)

        # add a fourth and fifth column as expected in writeDot (bug report 39)
        sifFile<-cbind(sifFile,BStimes)
        sifFile<-cbind(sifFile,weightsE)
        }
    else{
        #in this case there are AND gates and so we need to create dummy "and#" nodes
        sifFile<-matrix(0,nrow=4*length(reacOutput),ncol=5)
        nR<-1
        nANDs<-1

        for(i in 1:length(reacOutput)){
            if(length(reacInput[[i]]) == 1){
                sifFile[nR,1]<-reacInput[[i]]
                sifFile[nR,3]<-reacOutput[i]
                sifFile[nR,2]<-ifelse(any(modelComprExpanded$notMat[,i] == 1),-1,1)
                sifFile[nR,4]<-BStimes[i]
                sifFile[nR,5]<-weightsE[i]
                nR<-nR+1
                }
            else{
                for(inp in 1:length(reacInput[[i]])){
                    sifFile[nR,1]<-reacInput[[i]][inp]
                    sifFile[nR,3]<-paste("and",nANDs,sep="")
                    temp_indices = which(reacInput[[i]][inp]==rownames(modelComprExpanded$notMat))
                    sifFile[nR,2]<-ifelse(modelComprExpanded$notMat[temp_indices, i]==1,-1,1)
                    sifFile[nR,4]<-BStimes[i]
                    sifFile[nR,5]<-weightsE[i]
                    nR<-nR+1
                }
                sifFile[nR,1]<-paste("and",nANDs,sep="")
                sifFile[nR,3]<-reacOutput[i]
                sifFile[nR,2]<-1
                sifFile[nR,4]<-BStimes[i]
                sifFile[nR,5]<-weightsE[i]
                nANDs<-nANDs+1
                nR<-nR+1
                }
            }
        sifFile<-sifFile[1:(nR-1),]
        EApresent<-apply(sifFile[,1:3],1,createReac)
        EAweights<-cbind(EApresent,sifFile[,5])
        EApresent<-cbind(EApresent,sifFile[,4])
    }
    #this mini function makes edge attributes in the format e1 (sign) e2 = attr

    makeEA<-function(x){
        ea<-paste(x[1],"=",x[2])
        return(ea)
        }
    EApresent<-apply(EApresent,1,makeEA)
    EAweights<-apply(EAweights,1,makeEA)

#this is the edge attribute matrix that contains, for each edge, whether it is
#absent from the model (0), present at t1(1) or present at t2(2)
    return(list(EApresent=EApresent,EAweights=EAweights,sifFile=sifFile))
    }

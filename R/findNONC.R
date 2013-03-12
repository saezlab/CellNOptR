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
# $Id: findNONC.R 3155 2013-01-09 15:24:58Z cokelaer $
findNONC<-function(model,indexes,verbose=FALSE){

    #use the floyd warshall algorithm implemented in the package RBGL
    #eg floyd.warshall.all.pairs.sp(coex) that returns a distance matrix
    #for that I need to create a graphNEL object
    modeltoGraphNEL<-function(model){

        vertices<-model$namesSpecies
        edgeList<-vector("list", length=length(vertices))
        names(edgeList)<-vertices
        #edgeList is a list with an element edges and an element weights, both of the same size
        #where each element    in edgeL refers to one node of the graph and contains a vector made of
        #the indexes of the nodes to which the node in question has an edge

        for(sp in 1:length(vertices)){
            reacs<-which(model$interMat[sp,] == -1)

            if(length(reacs) == 0){
                targets<-sp

                }else{
                    targets<-which(model$interMat[,reacs[1]] == 1)

                    if(length(reacs) > 1){

                        for(r in 2:length(reacs)){
                            targets<-c(targets, which(model$interMat[,reacs[r]] == 1))
                            }

                        }
                    }
    #If the node doesn't have a target, I set it a dummy interaction with itself because
    #otherwise the graphNEL object can't be created. This won't change whether there exists
    #a path to it from a signal/cue or not
            edgeList[[sp]] <- list(edges=targets, weights=rep(1,length(targets)))
            }

        graph<-new("graphNEL", nodes=vertices, edgeL=edgeList,edgemode="directed")
        return(graph)
        }

    #Now I use the function above to create a GraphNEL object
    graphModel<-modeltoGraphNEL(model)

    #now I can run the floyd warshall algo
    distMatrix<-floyd.warshall.all.pairs.sp(graphModel)

    #Find the species that are non-observable, i.e. there are no path from that species to any
    #of the signals
    #And the species that are not controllable, ie there are no paths from any of the cues
    #to those species
    ncno<-0

    for(sp in 1:length(model$namesSpecies)){
    #if all the distances from this species to signals are infinite and the species is not a signal
    #then this species is not observable
        if(
            all(is.infinite(distMatrix[sp,indexes$signals]))  &&
            !any(indexes$signals == sp)){

            ncno<-c(ncno,sp)

            }

        if(
            all(is.infinite(distMatrix[c(indexes$stimulated,indexes$inhibited),sp])) &&
            !any(c(indexes$stimulated,indexes$inhibited) == sp)){

            ncno<-c(ncno,sp)

            }
        }
    ncno<-ncno[-1]
    ncno<-unique(ncno)

    #check that none of the cues or signals are in there
    cs<-intersect(unlist(indexes),ncno)

    if(length(cs != 0)){
        ncno<-ncno[-match(cs,ncno)]
        }

    if(verbose){
        print(paste(
            "The following species are not observable and/or not controllable:",
            toString(model$namesSpecies[ncno])))
        }
    return(ncno)
}


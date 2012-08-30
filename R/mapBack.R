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
#  CNO website: http://www.cellnopt.org
#
##############################################################################
# $Id$
mapBack <-
  function(model, PKN, bString){
    #the mapback for each link A->B in the compressed model is done looking
    #at the PKN as a graph and considering a subgraph of it including only
    #node A, node B and compressed nodes. All paths going from A to B in this
    #subnetwork are marked as active in the PKN.

    bStringPKN<-rep(0, length(PKN$reacID))

    #consider the PKN as a graph
    PKNgraph<-sif2graph(model2sif(model=PKN))

    #want to map back only selected links (1 in the bString)
    links2map<-model$reacID[as.logical(bString)]

    #function to search all paths connecting 2 nodes (from StartNode to EndNode)
    #in a given graph (gg)
    SearchPath<-function(gg,StartNode,EndNode){

      nod<-nodes(gg)
      edg<-edges(gg)

      #PathLis: list where all Paths starting from StartNode will e saved
      PathList<-list()
      ix<-1
      PathList[[ix]]<-StartNode

      #the search is implemented using a queue
      #initialization of the queue with all nodes reachale from the StartNode
      queue<-edg[[which(nod==StartNode)]]
      AllQueue<-queue  #this is just a backup al the whole queue for check
      #new visited paths are added
      for (k in 1:length(queue)){
        PathList[[ix+k]]<-c(PathList[[ix]],queue[k])
      }
      ix<-ix+1 #this index points always to the path walked to reach the current first element of the queue

      # remove nodes already visited to avoid problems with cyclic graphs
      gg<-removeNode(StartNode,gg)
      nod<-nodes(gg)
      edg<-edges(gg)

      # untill the queue is not empty
      while (length(queue)>0){
        if (any(which(nod==queue[1]))){
          #save in tmp nodes that can be reached from the first element of the queue
          tmp<-edg[[which(nod==queue[1])]]

          #update the queue with the new nodes at the end and save the corresponding paths walked to reach them
          if (length(tmp)>0){
            queue<-c(queue, tmp)
            AllQueue<-c(AllQueue, tmp)
            for (k in 1:length(tmp)){
              PathList[[length(PathList)+1]]<-as.vector(c(PathList[[ix]],tmp[k]))
            }
          }

          # already visited nodes are removed from the graph to avoid problems with cyclic graphs
          gg<-removeNode(queue[1],gg)
          nod<-nodes(gg)
          edg<-edges(gg)
        }
        # and the first element is removed from the queue
        queue<-queue[-1]
        ix<-ix+1
      }
      PathOK<-PathList[which( sapply(PathList,tail,1)==EndNode)]
      return(PathOK)
    }

    #for each hyperedge in the optimised model...
    for (i in 1:length(links2map)){
      x<-links2map[i]
      x<-gsub("!", "", x)
      x<-unlist(strsplit(x, "=", fixed=TRUE))
      EndNode<-x[2]
      # can contain more than 1 element in case of AND nodes
      StartNodes<-unlist(strsplit(x[1], "+", fixed=TRUE))
      for (j in 1:length(StartNodes)){
        StartNode<-StartNodes[j]

        #consider a subgraph including only StartNode, EndNode and compressed nodes
        okNodes<-c(model$speciesCompressed, EndNode, StartNode)
        noNodes<-setdiff(nodes(PKNgraph), okNodes)

        #this is the subgraph
        gg<-removeNode(noNodes, PKNgraph)

        #search this subgraph for paths between StartNode and EndNode
        PathOK<-SearchPath(gg=gg,StartNode=StartNode,EndNode=EndNode)

        for (j1 in 1:length(PathOK)){
          for (j2 in 2:length(PathOK[[j1]])){
            node1<-PathOK[[j1]][j2-1]
            node2<-PathOK[[j1]][j2]
            reacPKN<-paste(node1,node2,sep="=")
            bStringPKN[grep(reacPKN,PKN$reacID)]<-1
          }
        }

      }
    }


    findOutput<-function(x){
      sp<-which(x == 1)
      sp<-PKN$namesSpecies[sp]
    }
    reacOutput<-apply(PKN$interMat,2,findOutput)
    findInput<-function(x){
      sp<-which(x == -1)
      sp<-PKN$namesSpecies[sp]
    }
    reacInput<-apply(PKN$interMat,2,findInput)

    isNeg<-function(x){
      isNegI<-any(x == 1)
      return(isNegI)
    }
    inpSign<-apply(PKN$notMat,2,isNeg)
    inpSign<-!inpSign
    inpSign[inpSign]<-1
    inpSign[!inpSign]<--1

    sifFile<-cbind(reacInput,inpSign,reacOutput)
    createReac<-function(x){
      r<-paste(x[1]," (",x[2],") ",x[3],sep="")
      return(r)
    }
    EApresent<-apply(sifFile,1,createReac)
    EApresent<-cbind(EApresent,bStringPKN)

    makeEA<-function(x){
      ea<-paste(x[1],"=",x[2])
      return(ea)
    }
    EApresent<-apply(EApresent,1,makeEA)

    write.table(
      EApresent,
      file="weightsPKN.EA",
      row.names=FALSE,col.names="Weights",quote=FALSE,sep="\t")

    return(bStringPKN)
  }

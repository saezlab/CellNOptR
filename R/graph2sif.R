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
graph2sif<-function(graph,writeSif=FALSE, filename="Graph"){
  
  ed<-edges(graph)
  no<-nodes(graph)
  sifFile<-matrix(ncol=3)
  
  for (i in 1:length(no)){
    node1<-no[i]
    if (length(ed[[i]])>0){
      for (j in 1:length(ed[[i]])){
        node2<-ed[[i]][j]
        sig<-sign(as.numeric(edgeData(graph,no[i],ed[[i]][j], "weight")))
        sifFile<-rbind(sifFile,c(node1,sig,node2))
      }
    }
  }
  
  sifFile<-sifFile[-1,]
  
  if (writeSif==TRUE){
    filename<-paste(filename, ".sif", sep="")
    write.table(sifFile, file=filename, row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
  }

  return(sifFile)
}

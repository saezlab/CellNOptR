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
sif2graph<-function(sif){
	
  #if the input is a character it shoud be the name ot the sif file
  #otherwise a matrix in the sif format
  if (is.vector(sif) && (typeof(sif) == "character")){
	sif = read.table(sif) 
  }
	
  # build the unique vertices from the column 1 and 3 of the SIF file
  vertices = unique(c(as.character(sif[,1]), as.character(sif[,3])))
  # some aliases
  v1 = sif[,1]
  v2 = sif[,3]
  edges = as.numeric(sif[,2])

  l = length(vertices) - 1
  g <- new("graphNEL", nodes=vertices, edgemode="directed")
  #weights = rep(1, l)
  weights = edges
  for (i in 1:length(v1)){
    g <- addEdge(as.character(v1[i]), as.character(v2[i]), g, weights[i])
  }
  return(g)
  
}

#
#  This file is part of the CNO software
#
#  Copyright (c) 2019 - SaezLab - Heidelberg Universit
#
#  File author(s): CNO developers
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  CNO website: http://www.cellnopt.org
#
##############################################################################

link_midas_colnames_to_namesSpecies <- function(model, reducedMidas){
  
  linker = c()
  if (!is.null(colnames(reducedMidas))) { 
    for(i in 1:length(colnames(reducedMidas))){
      if(str_sub(colnames(reducedMidas)[i],-1)=="i"){
        linker[i] <- which(model$namesSpecies==
                            substr(strsplit(
                              colnames(reducedMidas)[i],":")[[1]][2],
                              1,nchar(strsplit(colnames(reducedMidas)[i],
                                               ":")[[1]][2])-1))
      } else{  ##for everything else
        linker[i] <- which(model$namesSpecies==
                            strsplit(colnames(reducedMidas)[i],":")[[1]][2])
      }
    }
  } else {
    linker <- 1:length(model$namesSpecies)
  }
  return(linker)
  
}
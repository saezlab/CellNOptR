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

createILPBitstringAll<- function(cplexSolutionFileName,
                                 y_vector,
                                 binary_variables){
  
  #
  bitStringAll <- list()
  
  #
  readOutBinariesAll<- function(cplexSolutionFileName, binary_variables){
    
    cplexSolutionData <- xmlParse(cplexSolutionFileName)
    cplexSolution <- xmlToList(cplexSolutionData, simplify = TRUE)
    binariesAll <- list()
    for(i in 1:length(cplexSolution)){
      
      if(all(names(cplexSolution[[i]])=="variable")){
        
        bitstring <- c()
        binaries <- c()
        for(j in 1:length(cplexSolution[[i]])){
          
          binaries[which(binary_variables[[1]]==cplexSolution[[i]][[j]][1])] = 
            cplexSolution[[i]][[j]][3]
          
        }
        
        binariesAll[[length(binariesAll)+1]] <- binaries
        
      }
      
    }
    
    return(binariesAll)
    
  }
  
  #
  binaries <- readOutBinariesAll(cplexSolutionFileName, binary_variables)
  for(i in 1:length(binaries)){
    bitString <- c()
    bitString <- binaries[[i]][1:length(y_vector)]
    bitStringAll[[length(bitStringAll)+1]] <- bitString
  }
  
  #
  bitStringILPAll <- bitStringAll
  for(i in 1:length(bitStringILPAll)){
    
    bitStringILPAll[[i]] <- round(as.numeric(bitStringILPAll[[i]]))
    
  }
  bitStringILPAll <- unique(bitStringILPAll)
  idx <- c()
  for(i in 1:length(bitStringILPAll)){
    
    if(sum(as.numeric(bitStringILPAll[[i]]))==0){
      
      idx <- c(idx, i)
      
    }
    
  }
  if(length(idx) > 0){
    
    bitStringILPAll <- bitStringILPAll[-idx]
    
  }
  
  return(bitStringILPAll)
  
}
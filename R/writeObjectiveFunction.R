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

writeObjectiveFunction <- function(model, 
                                   midasExperimentPart, 
                                   y_vector=y_vector,
                                   accountForModelSize = TRUE, 
                                   sizeFac = 0.000001, 
                                   meansOfMeasurements_at_t0, 
                                   method = "quadratic" ){
  
  # creating auxilliary objects
  reducedMidas <- midasExperimentPart
  whichObjective <- method
  objective1 <- c()
  
  # writing the objective function
  linker_vector <- link_midas_colnames_to_namesSpecies(model, reducedMidas)
  weight_of_unstimulated_experiments <- if(!is.null(dim(reducedMidas)[1])){dim(
    reducedMidas)[1]+1}else{1}
  constant <- 0
  factors <- "obj:\t"
  if(whichObjective == "linear"){
    for(i in 1:dim(reducedMidas)[1]){
      
      for(j in 1:dim(reducedMidas)[2]){
        if(!is.na(reducedMidas[i,j])){ 
          ##term1:
          if(i ==1){
            term1 <- weight_of_unstimulated_experiments * 
              meansOfMeasurements_at_t0[j]
          }else{
            term1 <- reducedMidas[i,j]
          }
          ##term2:
          if(i ==1){
            term2 <- weight_of_unstimulated_experiments * 
              (1-2*meansOfMeasurements_at_t0[j])
          }else{
            term2 <- 1-2*reducedMidas[i,j]
          }
          
          binary_index <- length(model$reacID) + 
            (i-1)*length(model$namesSpecies) + 
            linker_vector[j]  
          term3 <- paste0("xb", binary_index)
          constant <- constant + term1
          
          if(term2<0){
            factors  <- paste(factors,term2,term3)
          }else{
            factors <- paste(factors, "+", term2, term3)
          }
          
        }
        
        
      }
      if(accountForModelSize){
        for(m in 1:length(y_vector)){
          factors <- paste0(factors, " + ", sizeFac, " xb", m )
        }
      }
      
    }
    # here, the model size is accounted for:
    for(m in 1:length(y_vector)){
      factors <- paste0(factors, " + ", sizeFac, " xb", m )
    }
    
    out <- list(constant, factors) #see comment in the functon definition above
    
    
  } else if(whichObjective =="quadratic")
  {
    for(i in 1:dim(reducedMidas)[1]){
      
      for(j in 1:dim(reducedMidas)[2]){
        
        if(!is.na(reducedMidas[i, j])){
          
          binary_index <- length(model$reacID) + 
            (i-1)*length(model$namesSpecies) + linker_vector[j]
          #real term1:
          term1 <- paste0("xb", binary_index)
          #term2:
          if(reducedMidas[i, j]>=0){
            term2 <- paste0(2*abs(reducedMidas[i,j]), " xb", binary_index)
          } else {
            term2 <- paste0(0, " xb", binary_index)
          }
          #real term3:
          term3 <- reducedMidas[i,j] *reducedMidas[i,j]
          
          factors  <- paste(factors,"+", term1,"-", term2)#,"+", term3)
          
        }
        
        
      }
      
      if(accountForModelSize){
        for(m in 1:length(y_vector)){
          factors <- paste0(factors, " + ", sizeFac, " xb", m )
        }
      }
      
      out <- list(constant, factors)
    }
  }
  
  return(out)
  
}
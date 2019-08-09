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
                                   sizeFac = .000001, 
                                   meansOfMeasurements_at_t0, 
                                   method = "quadratic" ){ # function returns 2 values in a list: 1.) part of the objetive that only consists of absolute values, 2.) constant offset.
  
  #
  reducedMidas = midasExperimentPart
  whichObjective = method
  objective1 = c()
  
  # need link between column header and species in the binary_variables vector
  link_midas_colnames_to_namesSpecies <- function(model, reducedMidas){
    
    linker = c()
    if (!is.null(colnames(reducedMidas))) { # if there is no column name (single experiment) - added 13.02.18
      for(i in 1:length(colnames(reducedMidas))){
        if(str_sub(colnames(reducedMidas)[i],-1)=="i"){# for the inhibitors in the treatement section of the  midas file
          linker[i] = which(model$namesSpecies==substr(strsplit(colnames(reducedMidas)[i],":")[[1]][2],1,nchar(strsplit(colnames(reducedMidas)[i],":")[[1]][2])-1))
        } else{  #for everything else
          linker[i] = which(model$namesSpecies==strsplit(colnames(reducedMidas)[i],":")[[1]][2])
        }
      }
    } else {
      linker <- 1:length(model$namesSpecies)
    }
    return(linker)
    
  }
  
  linker_vector = link_midas_colnames_to_namesSpecies(model, reducedMidas)
  weight_of_unstimulated_experiments = if(!is.null(dim(reducedMidas)[1])){dim(reducedMidas)[1]+1}else{1} # this is N+1 / if clause added 13.02.18 
  constant = 0
  factors = "obj:\t"
  if(whichObjective == "linear"){
    for(i in 1:dim(reducedMidas)[1]){
      
      for(j in 1:dim(reducedMidas)[2]){
        if(!is.na(reducedMidas[i,j])){ #if measurements are NA (= not measured) they cannot appear in the objective function
          #term1:
          if(i ==1){
            term1 = weight_of_unstimulated_experiments * meansOfMeasurements_at_t0[j]
          }else{
            term1 = reducedMidas[i,j]
          }
          #term2:
          if(i ==1){
            term2 = weight_of_unstimulated_experiments * (1-2*meansOfMeasurements_at_t0[j])
          }else{
            term2 = 1-2*reducedMidas[i,j]
          }
          #term3:
          #here, the number of the binary is needed - procedure: use linker_vector to check, which number the respective species has in the namesSpecies Vector, then add 
          #(1.) the number of reactions (=offset), 
          #(2.) the number of experiments that have been evaluated before 
          #(3.) the index of the species evaluated in this experiment, 
          #(4.) (-1) since the binaries start with index 0 instead of 1.
          #dim(reducedMidasFile)[1]
          binary_index = length(model$reacID) + (i-1)*length(model$namesSpecies) + linker_vector[j]  # -1 is used since the binaries in Mitsos code start with xb0 - can be taken out after prove of principle
          term3 = paste0("xb", binary_index)#, "    ", binary_variables[[2]][binary_index])
          constant = constant + term1
          
          #print(paste0("Variable: xb", binary_index, " term1: ", term1, " term2: ", term2, " term3: ", term3, " meansOfMeasurements_at_t0: ", meansOfMeasurements_at_t0[j]))
          
          if(term2<0){
            factors  = paste(factors,term2,term3)
          }else{
            factors = paste(factors, "+", term2, term3)
          }
          
        }
        
        
      }
      if(accountForModelSize){
        for(m in 1:length(y_vector)){
          factors = paste0(factors, " + ", sizeFac, " xb", m )
        }
      }
      
    }
    # here, the model size is accounted for:
    for(m in 1:length(y_vector)){
      factors = paste0(factors, " + ", sizeFac, " xb", m )
    }
    
    out = list(constant, factors) #see comment in the functon definition above
    
    
  } else if(whichObjective =="quadratic")
  {
    for(i in 1:dim(reducedMidas)[1]){
      # for boolean predictions, the square difference between measurement and prediction can be calculated as
      # (x-x_m)^2 = x -   2*x*x_m +   x_m*x_m
      #           term_1    term2       term3
      for(j in 1:dim(reducedMidas)[2]){
        
        if(!is.na(reducedMidas[i, j])){
          
          binary_index = length(model$reacID) + (i-1)*length(model$namesSpecies) + linker_vector[j]  # -1 is used since the binaries in Mitsos code start with xb0 - can be taken out after prove of principle
          #real term1:
          term1 = paste0("xb", binary_index)
          #term2:
          if(reducedMidas[i, j]>=0){
            term2 = paste0(2*abs(reducedMidas[i,j]), " xb", binary_index)
          } else {
            term2 = paste0(0, " xb", binary_index)
          }
          #real term3:
          term3 = reducedMidas[i,j] *reducedMidas[i,j]
          
          factors  = paste(factors,"+", term1,"-", term2)#,"+", term3)
          
        }
        
        
      }
      
      
      out = list(constant, factors) #see comment in the functon definition above
    }
  }
  
  return(out)
  
}

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

write_bounds <- function(model,
                         midasTreatmentPart,
                         y_vector,
                         binary_variables){
                             
  #package stringr is needed in order to run this function. otherwise problems in the str_sub() command
  # this function returns the bounds vector which contains the boundaries for each binary varible.
  # procedure: create a vector that has the default binary bounds 0 <= xb_i <= 1
  # then, overwrite the elements of the vector that are fixed to a number (e.g. Treatment 1/0 or inhibitor 1) with the respective value.
  
  #
  writeBoundsForEdges_y_i<- function(y_vector){
    y_bound_vector = c()
    for(i in 1:length(y_vector)){
      y_bound_vector[i] = paste0("0 <= xb", i, " <= 1")
    }
    return(y_bound_vector)
  }
  
  #
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
  
  y_bounds_vector = writeBoundsForEdges_y_i(y_vector)
  linkerVector = link_midas_colnames_to_namesSpecies(model, midasTreatmentPart)
  bounds = c(paste0("0<= xb",1:(length(binary_variables[[1]])), " <= 1"))#, "   ", binary_variables[[2]][1:length(binary_variables[[3]])]))
  treatmentMatrix <- midasTreatmentPart
  for(i in 1:dim(midasTreatmentPart)[1]){
    for(j in 1:dim(midasTreatmentPart)[2]){
      if(midasTreatmentPart[i,j] == 1 && !(str_sub(colnames(midasTreatmentPart)[j], -1) == "i")){
        bounds[length(y_vector) + (i-1)*length(model$namesSpecies)+linkerVector[j]] = paste0("xb", length(y_vector) + (i-1)*length(model$namesSpecies)+linkerVector[j], " = 1")
      } else if(midasTreatmentPart[i,j] == 0 && !(str_sub(colnames(midasTreatmentPart)[j], -1) == "i")){
        bounds[length(y_vector) + (i-1)*length(model$namesSpecies)+linkerVector[j]] = paste0("xb", length(y_vector) + (i-1)*length(model$namesSpecies)+linkerVector[j], " = 0")
        
      }
      else if(midasTreatmentPart[i,j] == 1 && str_sub(colnames(midasTreatmentPart)[j], -1) == "i"){
        bounds[length(y_vector) + (i-1)*length(model$namesSpecies)+linkerVector[j]] = paste0("xb", length(y_vector) + (i-1)*length(model$namesSpecies)+linkerVector[j], " = 0")
      } else{
        #nothing. stays at it was.
        #bounds[length(y_vector) + (i-1)*length(model$namesSpecies) + linkerVector[j]-1] =
      }
    }
  }
  
  return(bounds)
}

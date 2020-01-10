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
  
  # writing of the boundaries
  y_bounds_vector = writeBoundsForEdges_y_i(y_vector)
  linkerVector = link_midas_colnames_to_namesSpecies(model, midasTreatmentPart)
  bounds = c(paste0("0<= xb",1:(length(binary_variables[[1]])), " <= 1"))
  treatmentMatrix <- midasTreatmentPart
  for(i in 1:dim(midasTreatmentPart)[1]){
    for(j in 1:dim(midasTreatmentPart)[2]){
      if(midasTreatmentPart[i,j] == 1 && 
         !(str_sub(colnames(midasTreatmentPart)[j], -1) == "i")){
        bounds[length(y_vector) + 
                 (i-1)*length(model$namesSpecies)+
                 linkerVector[j]] =
          paste0("xb", length(y_vector) + 
                   (i-1)*length(model$namesSpecies)+linkerVector[j], " = 1")
      } else if(midasTreatmentPart[i,j] == 0 && 
                !(str_sub(colnames(midasTreatmentPart)[j], -1) == "i")){
        bounds[length(y_vector) + 
                 (i-1)*length(model$namesSpecies)+
                 linkerVector[j]] = 
          paste0("xb", length(y_vector) + 
                   (i-1)*length(model$namesSpecies)+linkerVector[j], " = 0")
        
      }
      else if(midasTreatmentPart[i,j] == 1 && 
              str_sub(colnames(midasTreatmentPart)[j], -1) == "i"){
        bounds[length(y_vector) + 
                 (i-1)*length(model$namesSpecies)+
                 linkerVector[j]] = 
          paste0("xb", length(y_vector) + 
                   (i-1)*length(model$namesSpecies)+linkerVector[j], " = 0")
      } else{
        ## nothing. stays at it was.
      }
    }
  }
  
  return(bounds)
}
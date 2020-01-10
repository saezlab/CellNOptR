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

writeConstraint3 <- function(model, midasExperimentPart){
  
  reducedMidas = midasExperimentPart
  constraints = c()
  counter = 1
  for(i in 1:dim(reducedMidas)[1]){
    for(j in 1:length(model$reacID)){
      newConstraint = paste0("c", counter, ":\t", 
                             "xb", 
                             (length(model$reacID)+
                                dim(reducedMidas)[1]*
                                length(model$namesSpecies)+
                                length(model$reacID)*(i-1)+j), 
                             " - xb" , j, " <= 0")
      constraints[counter] =newConstraint
      counter= counter+1
    }
  }
  
  return(constraints)
  
}
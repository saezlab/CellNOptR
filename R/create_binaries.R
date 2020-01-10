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

create_binaries <- function(model, 
                            midas, 
                            numberOfExperiments, 
                            y_vector){
  
  ## creating auxilliary variables
  numbers <- c()
  variables <- c()
  identifiers <- c()
  
  ## Creating binary variables for the edges
  y_binaries <- create_binary_variables_for_y_vector(y_vector)
  numbers <- y_binaries[[1]]
  variables <- y_binaries[[2]]
  identifiers <- y_binaries[[3]]
  
  ## creating binary variables for the species
  for(i in 1:numberOfExperiments){
    binList <- 
      create_binary_variables_for_all_species_in_experiment_k_i(model, i)
    numbers= append(numbers, binList[[1]])
    variables <- append(variables, binList[[2]])
    identifiers <- append(identifiers, binList[[3]])
  }
  for(i in 1:numberOfExperiments){
    binList <- 
      create_binary_variables_for_all_reactions_z_i_in_experiment_k_i(model, 
                                                                      i)
    numbers= append(numbers, binList[[1]])
    variables <- append(variables, binList[[2]])
    identifiers <- append(identifiers, binList[[3]])
  }
  
  ## combining and returning all the binary variables
  binaries <- list(c(paste0("xb",1:length(variables))), variables, identifiers)
  return(binaries)
  
}
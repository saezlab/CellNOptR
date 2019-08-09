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
  
  #
  create_binary_variables_for_y_vector <- function(y_vector){
    
    variables = c()
    for(i in 1:length(y_vector)){
      variables[i] = paste0("y_",i)
    }
    binary_variables_list = list(1:(length(y_vector)),variables, paste("reaction",y_vector))
    return(binary_variables_list)
    
  }
  
  #
  create_binary_variables_for_all_species_in_experiment_k_i <- function(model, 
                                                                        numberOfExperiment){ # k means experiment index
    
    number_of_species = length(model$namesSpecies)
    identifier = c()
    variables = c()
    for(i in 1:number_of_species){
      identifier[i] = paste("species", model$namesSpecies[i], "in experiment k =", numberOfExperiment)
      variables[i] = paste0("x_",i,"^",numberOfExperiment)
    }
    binary_variables_list = list(0:length(model$namesSpecies),variables, identifier)
    return(binary_variables_list)
    
  }
  
  #
  create_binary_variables_for_all_reactions_z_i_in_experiment_k_i <- function(model, 
                                                                              numberOfExperiment){ # k means experiment index
    number_of_reactions = length(model$reacID)
    identifier = c()
    variables = c()
    for(i in 1:number_of_reactions){
      identifier[i] = paste("reaction", model$reacID[i], "in experiment k =", numberOfExperiment)
      variables[i] = paste0("z_",i,"^",numberOfExperiment)
    }
    binary_variables_list = list(0:length(model$reacID),variables, identifier)
    return(binary_variables_list)
  }
  
  #
  numbers = c()
  variables = c()
  identifiers = c()
  
  #
  y_binaries = create_binary_variables_for_y_vector(y_vector)
  numbers = y_binaries[[1]]
  variables = y_binaries[[2]]
  identifiers = y_binaries[[3]]
  
  #
  for(i in 1:numberOfExperiments){
    binList = create_binary_variables_for_all_species_in_experiment_k_i(model, i)
    numbers= append(numbers, binList[[1]])
    variables = append(variables, binList[[2]])
    identifiers = append(identifiers, binList[[3]])
  }
  for(i in 1:numberOfExperiments){
    binList = create_binary_variables_for_all_reactions_z_i_in_experiment_k_i(model, i)
    numbers= append(numbers, binList[[1]])
    variables = append(variables, binList[[2]])
    identifiers = append(identifiers, binList[[3]])
  }
  
  #
  binaries = list(c(paste0("xb",1:length(variables))), variables, identifiers)
  return(binaries)
  
}
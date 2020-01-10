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

create_binary_variables_for_all_species_in_experiment_k_i <- 
  function(model, numberOfExperiment){
    
    number_of_species <- length(model$namesSpecies)
    identifier <- c()
    variables <- c()
    for(i in 1:number_of_species){
      identifier[i] <- paste("species", model$namesSpecies[i], 
                            "in experiment k =", 
                            numberOfExperiment)
      variables[i] <- paste0("x_",i,"^",numberOfExperiment)
    }
    binary_variables_list <- list(0:length(model$namesSpecies),
                                 variables, 
                                 identifier)
    return(binary_variables_list)
    
}
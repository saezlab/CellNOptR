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

create_binary_variables_for_all_reactions_z_i_in_experiment_k_i <- 
  function(model, numberOfExperiment){ 
    number_of_reactions <- length(model$reacID)
    identifier <- c()
    variables <- c()
    for(i in 1:number_of_reactions){
      identifier[i] <- paste("reaction", model$reacID[i], "in experiment k =", 
                            numberOfExperiment)
      variables[i] <- paste0("z_",i,"^",numberOfExperiment)
    }
    binary_variables_list <- list(0:length(model$reacID),variables, identifier)
    return(binary_variables_list)
}
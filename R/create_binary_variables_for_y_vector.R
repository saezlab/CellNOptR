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

create_binary_variables_for_y_vector <- function(y_vector){
  
  variables = c()
  for(i in 1:length(y_vector)){
    variables[i] <- paste0("y_",i)
  }
  binary_variables_list <- list(1:(length(y_vector)),
                               variables, 
                               paste("reaction",y_vector))
  return(binary_variables_list)
  
}
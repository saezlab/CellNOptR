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

writeBoundsForEdges_y_i<- function(y_vector){
  y_bound_vector <- c()
  for(i in 1:length(y_vector)){
    y_bound_vector[i] <- paste0("0 <= xb", i, " <= 1")
  }
  return(y_bound_vector)
}
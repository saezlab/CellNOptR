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

writeFile <- function(objectiveFunction,
                      constraints,
                      bounds,
                      binaries){
  
  data = "testFile.lp"
  write("enter Problem", data)
  write("", data, append = TRUE)
  write("Minimize", data, append = TRUE)
  
  write(objectiveFunction[[2]], data, append = TRUE)
  
  write("Subject To", data, append = TRUE)
  write(constraints, data, append = TRUE)
  
  write("Bounds", data, append = TRUE)
  write(bounds, data, append = TRUE)
  
  write("Binaries",data, append = TRUE)
  write(binaries[[1]], data, append = TRUE)
  
  write("End", data, append = TRUE)
  write("optimize Problem", data, append = TRUE)
  write("write solution_cplex_primaryObjective.txt", data, append = TRUE)
  write("sol", data, append = TRUE)
  write("y", data, append = TRUE)
  write("quit", data, append = TRUE)
  
}
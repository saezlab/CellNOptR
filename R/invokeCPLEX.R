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

invokeCPLEX <- function(inputFileName, 
                        outputFileName, 
                        mipGap=mipGap, 
                        relGap = relGap, 
                        timelimit=timelimit, 
                        cplexPath = cplexPath, 
                        numSolutions = numSolutions, 
                        limitPop = limitPop, 
                        poolIntensity = poolIntensity, 
                        poolReplace = poolReplace){
  
  ## Writing and executing the cplex commands
  ## setwd(gsub(pattern = "cplex", replacement = "", x = cplexPath))
  writeCplexCommandFile(commandFileName = "cplexCommand.txt", inputFileName, 
                        outputFileName, mipGap=mipGap, relGap = relGap, 
                        timelimit=timelimit, cplexPath = cplexPath, 
                        numSolutions = numSolutions, limitPop = limitPop, 
                        poolIntensity = poolIntensity, 
                        poolReplace = poolReplace)
  if(file.exists(outputFileName)){
    file.remove(outputFileName)
  }
  # system(paste0(cplexPath, " -f cplexCommand.txt"))
  
  if (Sys.info()[1]=="Windows") {
    file.copy(from = cplexPath,to = getwd())
    system(paste0("cplex.exe -f cplexCommand.txt"))
    file.remove("cplex.exe")
  } else {
    system(paste0(cplexPath, " -f cplexCommand.txt"))
  }
  
}
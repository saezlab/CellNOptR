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
  
  # Function that writes the controll file for CPLEX. If needed, here additional options such as maximal computation time can be introduced
  writeCplexCommandFile <- function(commandFileName, 
                                    inputFileName, 
                                    outputFileName,
                                    mipGap=mipGap,
                                    relGap=relGap,
                                    timelimit=timelimit, 
                                    cplexPath=cplexPath, 
                                    numSolutions = numSolutions, 
                                    limitPop = limitPop, 
                                    poolIntensity = poolIntensity, 
                                    poolReplace = poolReplace){
    
    setwd(gsub(pattern = "cplex", replacement = "", x = cplexPath))
    data = commandFileName
    write(paste0("read ", inputFileName), data)
    write(paste0("set mip tolerances mipgap ",mipGap), data, append = TRUE)
    write(paste0("set mip pool relgap ",relGap), data, append = TRUE)
    write(paste0("set mip pool replace ", poolReplace), data, append = TRUE)
    write(paste0("set mip pool intensity ",poolIntensity), data, append = TRUE)
    write(paste0("set timelimit ",timelimit=timelimit), data, append = TRUE)
    write(paste0("set mip limits populate ", limitPop), data, append = TRUE)
    write(paste0("set mip pool capacity ", numSolutions), data, append = TRUE)
    write(paste0("populate"), data, append = TRUE)
    write(paste0("write ", outputFileName, " sol all"), data, append = TRUE)
    write("quit", data, append = TRUE)
    
  }
  
  setwd(gsub(pattern = "cplex", replacement = "", x = cplexPath))
  writeCplexCommandFile(commandFileName = "cplexCommand.txt", inputFileName, outputFileName, mipGap=mipGap, relGap = relGap, timelimit=timelimit, cplexPath = cplexPath, numSolutions = numSolutions, limitPop = limitPop, poolIntensity = poolIntensity, poolReplace = poolReplace)
  if(file.exists(outputFileName)){
    file.remove(outputFileName)
  }
  system(paste0(cplexPath, " -f cplexCommand.txt"))
  
}
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
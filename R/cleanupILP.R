#  Copyright (c) 2019 - SaezLab - Heidelberg University
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

cleanupILP <- function(){
  
  if(file.exists("testFile.lp")){
    file.remove("testFile.lp")
  }
  
  if(file.exists("cplexSolution.txt")){
    file.remove("cplexSolution.txt")
  }
  
  if(file.exists("cplex.log")){
    file.remove("cplex.log")
  }
  
  if(file.exists("cplexCommand.txt")){
    file.remove("cplexCommand.txt")
  }
  
  AllFiles <- list.files()
  CloneFiles <- which(grepl(pattern = "clone",x = AllFiles,fixed = TRUE))
  if (length(CloneFiles)>0) {
    for (counter in 1:length(CloneFiles)) {
      file.remove(AllFiles[CloneFiles[counter]])
    }
  }
  
}
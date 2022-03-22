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

ilpBinaryT2 <- function(cnolist,
                        model, 
                        sizeFac = 0.0001, 
                        mipGap = 0, 
                        relGap = 0, 
                        timelimit = 3600, 
                        cplexPath, 
                        method = "quadratic",
                        numSolutions = 100, 
                        limitPop = 500, 
                        poolIntensity = 0, 
                        poolReplace = 2){
  
  ## Initializing auxilliary objects for the analysis
  resILPAll <- list()
  exclusionList <- NULL
  cnolistReal <- cnolist
  tmpPath <- getwd()
  
  if(!file.exists(cplexPath)){
  	stop("User should provide a valid CPLEX path for using this function.")
  }
  
  ## Checking cnolist class
  if (is(CNOlist,"CNOlist")){
    writeMIDAS(CNOlist = cnolist, filename = "tempMD.csv", 
               timeIndices = c(1, 3), overwrite = TRUE)
    md <- readMIDAS(MIDASfile = "tempMD.csv")
    file.remove("tempMD.csv")
    cnolist <- compatCNOlist(object = cnolist)
  }
  
  ## Form of the objective function
  if(!(method%in%c("quadratic", "linear"))){
    print("'method' variable should be quadratic or linear. 
          Setting method = 'qadratic'")
    method <- "quadratic"
  }
  
  # Penalty factors as decimals
  options(scipen = 10)
  
  # Creating CNOlist object to use from the midas md table and checking if 
  ## signals in the cnolist and pkn match
  cnolist <- makeCNOlist(md,FALSE)
  
  startILP <- Sys.time()
  print("Creating LP file and running ILP...")
  resILP = createAndRunILP(model, md, cnolistReal, accountForModelSize = TRUE, 
                           sizeFac = sizeFac, mipGap=mipGap, relGap=relGap, 
                           timelimit=timelimit, cplexPath = cplexPath, 
                           method = method, numSolutions = numSolutions, 
                           limitPop = limitPop, poolIntensity = poolIntensity, 
                           poolReplace = poolReplace)
  endILP <- Sys.time()
  
  cleanupILP()
  
  return(resILP)
  
}

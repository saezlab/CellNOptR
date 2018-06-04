#
#  This file is part of the CNO software
#
#  Copyright (c) 2011-2018 - EBI/RWTH Aachen
#
#  File author(s): CNO developers (cno-dev@ebi.ac.uk)
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  CNO website: http://CellNOpt.org
#
##############################################################################
# $Id$

CNObootstrap <- function(model,CNOlistBS,BSmethod=1,N_bs=100,AddSD=0.05){
  
# Perform bootstrapping of dataset in CNOlist with permutation +/- with random sampling or resample from fixed mean and variance
# Input: CNOlist with the original dataset
# Output: CNOlist with the bootstrapped dataset
  
  CNOlistBSOrig <- CNOlistBS # Keep original
  
  # Generate new data to perform bootstrapping
  BSdataAll <- list()
  for (counter in 1:length(CNOlistBS$valueSignals)) {
    BSdataAll[[counter]] <- array(NA,dim = c(dim(CNOlistBS$valueSignals[[1]]),N_bs))
  }
  BSres <- list()
  
  rnorm2 <- function(n,mean,sd) { mean+sd*scale(rnorm(n)) }

  if (BSmethod == 1) { # Variant 1 - permutation/shuffling with/without replacement
    for (counter_BS in 1:N_bs) {
      for (counter in 1:length(CNOlistBS$valueSignals)) {
        for (counter2 in 1:ncol(CNOlistBS$valueSignals[[counter]])) {
          # Here shuffle the order of the measurement
          BSdataAll[[counter]][,,counter_BS][,counter2] <- CNOlistBS$valueSignals[[counter]][,counter2][sample(1:length(CNOlistBS$valueSignals[[counter]][,counter2]),size=length(CNOlistBS$valueSignals[[counter]][,counter2]),replace=TRUE)]      
        }
      }
    }
  } else if (BSmethod == 2) { # Variant 2 - sample from normal distribution
    if (is.null(CNOlistBS$valueVariances)) {
      CNOlistBS$valueVariances[[1]] <- CNOlistBS$valueSignals[[1]]
      CNOlistBS$valueVariances[[2]] <- matrix(0,nrow(CNOlistBS$valueVariances[[1]]),ncol(CNOlistBS$valueVariances[[1]]))
    }
    CNOlistBS$valueVariances[[2]] <- matrix(rep(AddSD,length(CNOlistBS$valueVariances[[2]])),nrow(CNOlistBS$valueVariances[[2]]),ncol(CNOlistBS$valueVariances[[2]]))
    # Assign function to perform sampling
    for (counter in 1:length(CNOlistBS$valueSignals)) {
      for (counter2 in 1:ncol(CNOlistBS$valueSignals[[counter]])) {
        for (counter3 in 1:nrow(CNOlistBS$valueSignals[[counter]])) {
          # Here re-sample from the normal distribution of mean and SD
          SampledData <- rnorm2(N_bs,CNOlistBS$valueSignals[[counter]][counter3,counter2],CNOlistBS$valueVariances[[counter]][counter3,counter2])
          for (counter4 in 1:length(SampledData)) {
            BSdataAll[[counter]][,,counter4][counter3,counter2] <- SampledData[counter4]
          }
        }
      }
    }
  }
  
  for (counter_BS in 1:N_bs) {
    
    print(paste0("Bootstrapping Round: ",counter_BS,"/",N_bs))
    
    CNOlistBS <- CNOlistBSOrig
    
    for (counter in 1:length(CNOlistBS$valueSignals)) {
      CNOlistBS$valueSignals[[counter]] <- BSdataAll[[counter]][,,counter_BS]
    }
    initBstring<-rep(1,length(model$reacID))
    ToyT1opt<-gaBinaryT1(CNOlist=CNOlistBS, model=model,initBstring=initBstring, verbose=FALSE)
    BSres[[counter_BS]] <- list("FitCost"=ToyT1opt$bScore,"FitParam"=ToyT1opt$bString)
    
  }
  
  # Process all fitting costs and fitted parameters
  AllFitCost <- NULL; AllFitParam <- NULL
  for (counter in 1:length(BSres)) {
    AllFitCost <- c(AllFitCost,BSres[[counter]]$FitCost)
    AllFitParam <- rbind(AllFitParam,BSres[[counter]]$FitParam)
  }
  
  pdf("FitCost_Bootstrapping.pdf")
  boxplot(x = AllFitCost, outpch = NA,main="Fitting Cost Bootstrapping") 
  stripchart(x = AllFitCost, 
             vertical = TRUE, method = "jitter", 
             pch = 21, col = "maroon", bg = "bisque", 
             add = TRUE) 
  dev.off()
  
  # Show figure on Knitr
  boxplot(x = AllFitCost, outpch = NA,main="Fitting Cost Bootstrapping") 
  stripchart(x = AllFitCost, 
             vertical = TRUE, method = "jitter", 
             pch = 21, col = "maroon", bg = "bisque", 
             add = TRUE) 
  
  pdf("FitParam_Bootstrapping.pdf")
  plotModel(model=model,CNOlist = CNOlistBSOrig,bString = colMeans(AllFitParam))
  dev.off()
  
  # Show figure on Knitr
  plotModel(model=model,CNOlist = CNOlistBSOrig,bString = colMeans(AllFitParam))
  
  return (list(AllFitCost=AllFitCost,AllFitParam=AllFitParam))
  
}
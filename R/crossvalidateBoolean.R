#' ## Crossvalidation study
#' crossvalidate
#' 
#' k-fold crossvalidation for Boolean model
#' 
#' Does a k-fold cross-validation for Boolean CellNOpt models. In k-iterations a 
#' fraction of the data is eliminated from the CNOlist. The model is trained 
#' on the 
#' remaining data and then the model predicts the held-out data. Then the 
#' prediction
#' accuracy is reported for each iteration. 
#' 
#' @param CNOlist Cnolist which contains all the experiments  
#' @param model a model prepared for the training  
#' @param nfolds number of folds - default is 10. Although nfolds can be as 
#' large as the sample size (leave-one-out CV),
#'  it is not recommended for large datasets.  
#' @param foldid an optional vector of values between `1` and `nfold` 
#' identifying what fold each observation is in. If supplied, `nfold` can be 
#' missing.  
#' @param type define the way to do the crossvalidation. The default is 
#' `type="datapoint"`, which assigns the data randomly into folds. 
#' The option `type="experiment"` uses whole experiments for crossvalidation 
#' (all data corresponding to a cue combination). The `type=observable` uses the
#' subset of nodes across all experiments for crossvalidation.  
#' @param ... further arguments are passed to gaBinaryT1  
#' @seealso \link{\code{gaBinaryT1}}  

crossvalidateBoolean <- function(CNOlist,model,nfolds=10,foldid=NULL, 
                                 type=c('datapoint','experiment','observable'),
                                 timeIndex = 2,parallel=FALSE, ...){
  
  if ((class(CNOlist)=="CNOlist")==FALSE){
    CNOlist = CellNOptR::CNOlist(CNOlist)
  }
  
  crossvalidate.call = match.call(expand.dots = TRUE)
  
  if(!(type%in%c('datapoint','experiment','observable'))){
    
    error("Wrong argument for type. It should either be 'datapoint', 
          'experiment' or 'observable'")
    
  } else {
    
    if(type=='datapoint'){
      N = prod(dim(CNOlist@signals[[timeIndex]]))
    }else if(type=='experiment'){
      N = nrow(CNOlist@signals[[timeIndex]])
    }else if(type=="observable"){
      N = ncol(CNOlist@signals[[timeIndex]])
    }
    
  }
  
  if (is.null(foldid)){
    
    foldid = sample(rep(seq(nfolds), length = N))
    
  } else {
    
    nfolds = max(foldid)
    
  }
  
  outlist = as.list(seq(nfolds))
  
  if (parallel) {
  	if (!requireNamespace("doParallel", quietly = TRUE)) {
  		stop("Package \"doParallel\" needed for this function to work. Please install it.",
  			 call. = FALSE)
  	}
  	
    outlist = foreach::foreach(i = seq(nfolds), .packages = c("CellNOptR")) %dopar%
    {
      whichI = foldid == i
      CNOlist.sub = CNOlist
      CNOlist.cv = CNOlist
      if(type=='datapoint'){
        CNOlist.sub@signals[[timeIndex]][whichI] = NA
        CNOlist.cv@signals[[timeIndex]][!whichI] = NA
      }else if(type=='experiment'){
        CNOlist.sub@signals[[timeIndex]][whichI,] = NA
        CNOlist.cv@signals[[timeIndex]][!whichI,] = NA
      }else if(type=="observable"){
        CNOlist.sub@signals[[timeIndex]][,whichI] = NA
        CNOlist.cv@signals[[timeIndex]][,!whichI] = NA
      }
      
      outlist = list()
      outlist$fit = gaBinaryT1(CNOlist = CNOlist.sub, 
                               model=model, timeIndex = timeIndex, verbose = FALSE)
      
      outlist$cvScore = computeScoreT1(CNOlist = CNOlist.cv,model = model,
                                       bString = outlist$fit$bString,
                                       timeIndex = timeIndex)
      outlist
    }
    
  }
  else {
    cat("fold \t fitScore \t cvScore \t nTolModels \t bString\n")
    
    for (i in seq(nfolds)) {
      whichI = foldid == i
      CNOlist.sub = CNOlist
      CNOlist.cv = CNOlist
      if(type=='datapoint'){
        CNOlist.sub@signals[[timeIndex]][whichI] = NA
        CNOlist.cv@signals[[timeIndex]][!whichI] = NA
      }else if(type=='experiment'){
        CNOlist.sub@signals[[timeIndex]][whichI,] = NA
        CNOlist.cv@signals[[timeIndex]][!whichI,] = NA
      }else if(type=="observable"){
        CNOlist.sub@signals[[timeIndex]][,whichI] = NA
        CNOlist.cv@signals[[timeIndex]][,!whichI] = NA
      }
      
      outlist[[i]] = list()
      outlist[[i]]$fit = gaBinaryT1(CNOlist = CNOlist.sub, 
                                    model=model, timeIndex = timeIndex, 
                                    verbose = FALSE, ...)
      
      outlist[[i]]$cvScore = computeScoreT1(CNOlist = CNOlist.cv,model = model,
                                            bString = outlist[[i]]$fit$bString,
                                            timeIndex = timeIndex)
      cat(i,"\t",outlist[[i]]$fit$bScore,"\t",outlist[[i]]$cvScore,"\t",
          length(outlist[[i]]$fit$stringsTolScores),"\t\t",
          paste(outlist[[i]]$fit$bString,collapse=","),"\n")
    }
  }
  
  cvScores = data.frame(folds=seq(nfolds), cvScore= 
                          unlist(lapply(outlist,function(x)x$cvScore)))
  fitScores = data.frame(folds=seq(nfolds), bScore= 
                           unlist(lapply(outlist,
                                         function(x)x$fit[c("bScore")])))
  bStrings = data.frame(folds=seq(nfolds), 
                        do.call("rbind",lapply(outlist,
                                               function(x)matrix(x$fit$bString,
                                                                 nrow=1))))
  colnames(bStrings) = c('folds', model$reacID)
  
  
  out = list(cvScores=cvScores,
             fitScores=fitScores,
             bStrings=bStrings,
             crossvalidate.call=crossvalidate.call,
             foldid=foldid)
  return(out)
}

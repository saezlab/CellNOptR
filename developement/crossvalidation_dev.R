#' # Crossvalidation for CellNOptR Boolean models  
#' author: Attila Gabor  
#' date: 25.05.2018  

library(CellNOptR)

#' ## MacNamara case study:
pknmodel = readSIF("caseStudies/PKN-ToyPB.sif")
cnodata = CNOlist("caseStudies/MD-ToyPB.csv")

#' original and preprocessed network 
plotModel(pknmodel,cnodata)
model = preprocessing(data = cnodata,model = pknmodel,compression = T,expansion = T)
plotModel(model,cnodata)

#' original CNOlist contains many timepoints, we use only a subset
plot(cnodata)
selectedTime = c(0,10)
cnodata_prep = cutCNOlist(cnodata,model = model,cutTimeIndices = which(!cnodata@timepoints %in% selectedTime))
plot(cnodata_prep)

#' optimise and show results
opt = gaBinaryT1(CNOlist = cnodata_prep,model = model,verbose = F)
print(opt$bString)
print(opt$bScore)

S = cutAndPlot(CNOlist = cnodata_prep,model = model,bStrings =list(opt$bString))

plotModel(model = model,CNOlist = cnodata_prep,bString = opt$bString)

#' ## Crossvalidation study
#' crossvalidate
#' 
#' k-fold crossvalidation for Boolean model
#' 
#' Does a k-fold cross-validation for Boolean CellNOpt models. In k-iterations a 
#' fraction of the data is eliminated from the CNOlist. The model is trained on the 
#' remaining data and then the model predicts the held-out data. Then the prediction
#' accuracy is reported for each iteration. 
#' 
#' @param CNOlist Cnolist which contains all the experiments  
#' @param model a model prepared for the training  
#' @param nfolds number of folds - default is 10. Although nfolds can be as large as the sample size (leave-one-out CV),
#'  it is not recommended for large datasets.  
#' @param foldid an optional vector of values between `1` and `nfold` identifying what fold each observation is in. If supplied, `nfold` can be missing.  
#' @param type define the way to do the crossvalidation. The default is 
#' `type="datapoint"`, which assigns the data randomly into folds. 
#' The option `type="experiment"` uses whole experiments for crossvalidation 
#' (all data corresponding to a cue combination). The `type=observable` uses the
#' subset of nodes across all experiments for crossvalidation.  
#' @param ... further arguments are passed to gaBinaryT1  
#' @seealso \link{\code{gaBinaryT1}}  

crossvalidate = function(CNOlist,model,nfolds=10,foldid, type=c('datapoint','experiment','observable'),timeIndex = 2,parallel=FALSE, ...){
	
	if ((class(CNOlist)=="CNOlist")==FALSE){
		CNOlist = CellNOptR::CNOlist(CNOlist)
	}
	
	crossvalidate.call = match.call(expand.dots = TRUE)
	
	type = match.arg(type)
	if(type=='datapoint'){
		N = prod(dim(CNOlist@signals[[timeIndex]]))
	}else if(type=='experiment'){
		N = nrow(CNOlist@signals[[timeIndex]])
	}else if(type=="observable"){
		N = ncol(CNOlist@signals[[timeIndex]])
	}
	
	if (missing(foldid)) 
		foldid = sample(rep(seq(nfolds), length = N))
	else nfolds = max(foldid)
	
	outlist = as.list(seq(nfolds))
	
	if (parallel) {
		require(doParallel)
		outlist = foreach(i = seq(nfolds), .packages = c("CellNOptR")) %dopar% 
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
										  model=model, timeIndex = timeIndex, verbose = F)
			
			
			# simRes = simulateTN(CNOlist = CNOlist.cv,model = model,bStrings = list(outlist[[i]]$fit$bString))
			
			outlist$cvScore = computeScoreT1(CNOlist = CNOlist.cv,model = model,bString = outlist$fit$bString,timeIndex = timeIndex)
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
									  model=model, timeIndex = timeIndex, verbose = F, ...)
			
			
			# simRes = simulateTN(CNOlist = CNOlist.cv,model = model,bStrings = list(outlist[[i]]$fit$bString))
			
			outlist[[i]]$cvScore = computeScoreT1(CNOlist = CNOlist.cv,model = model,bString = outlist[[i]]$fit$bString,timeIndex = timeIndex)
			# Score = getFit(simResults = simRes,
			# 	   CNOlist = CNOlist.cv,
			# 	   model = model,
			# 	   timePoint = 2,
			# 	   indexList = indexFinder(CNOlist.cv,model),
			# 	   nInTot=length(which(model$interMat == -1))
			# 	   )
			# nDataP <- sum(!is.na(CNOlist.cv@signals[[timeIndex]]))
			# outlist[[i]]$cvScore <- Score/nDataP
			# 
			cat(i,"\t",outlist[[i]]$fit$bScore,"\t",outlist[[i]]$cvScore,"\t",
				length(outlist[[i]]$fit$stringsTolScores),"\t\t",
				paste(outlist[[i]]$fit$bString,collapse=","),"\n")
		}
	}
	
	cvScores = data.frame(folds=seq(nfolds), cvScore= unlist(lapply(outlist,function(x)x$cvScore)))
	fitScores = data.frame(folds=seq(nfolds), bScore= unlist(lapply(outlist,function(x)x$fit[c("bScore")])))
	bStrings = data.frame(folds=seq(nfolds), do.call("rbind",lapply(outlist,function(x)matrix(x$fit$bString,nrow=1))))
	colnames(bStrings) = c('folds', model$reacID)
	
	
	out = list(cvScores=cvScores,
			   fitScores=fitScores,
			   bStrings=bStrings,
			   crossvalidate.call=crossvalidate.call,
			   foldid=foldid)
	return(out)
}


#' ## 10-fold crossvalidation using T1 data
#' We use only T1 data for crossvalidation, because data in the T0 matrix is not independent.
#' All rows of data in T0 describes the basal condition.

#' Crossvalidation produce some text in the command window:  
#' *Shall we keep it?*
system.time({R1 = crossvalidate(cnodata_prep, model=model, type="datapoint", nfolds=10)})[[3]]

library(doParallel)
registerDoParallel(cores=3)
system.time({R2=crossvalidate(cnodata_prep, model=model, type="datapoint", nfolds=10, parallel = TRUE)})[[3]]

#' we pick the fold that gave the worse CV score to see what the code is doing
# select fold
worseFold = which.max(R1$cvScores$cvScore)
whichI = R1$foldid == worseFold
# prepare training and CV data
timeIndex = 2
CNOlist.sub = cnodata_prep
CNOlist.cv = cnodata_prep

CNOlist.sub@signals[[timeIndex]][whichI] = NA
CNOlist.cv@signals[[timeIndex]][!whichI] = NA

#' Data used for training:
plot(CNOlist.sub)

#' Cross-validation dataset:
plot(CNOlist.cv)

#' We can see that the training data and the crossvalidation data overlaps regarding time point 0, but 
#' they are totally distinct at time point 1.  

# Training the model on the fitting data
fit = gaBinaryT1(CNOlist = CNOlist.sub, 
				  model=model, timeIndex = timeIndex, verbose = F)

print(fit$bScore)
print(fit$bString)

#' comparing to the bitString of the optimal model
print(opt$bString)

# if the bitStrings are different, plot the network for comparison. 
if(any(opt$bString!=fit$bString)){
	plotModel(model = model,CNOlist = CNOlist.sub,bString = fit$bString)
}

#' optimal model's score on the training set:
optScore = computeScoreT1(CNOlist = CNOlist.sub,model = model,bString = opt$bString)
print(optScore)

#' show the model fit on the training set:
S1 = cutAndPlot(CNOlist.sub, model, bStrings=list(fit$bString))


#' model predictions on the cross-validation set
S1 = cutAndPlot(CNOlist.cv, model, bStrings=list(fit$bString))


#' compute score for the cross-validation
cvScore = computeScoreT1(CNOlist = CNOlist.cv,model = model,bString = fit$bString)

#' score for cross-validation:
print(cvScore)




#' ## 10-fold crossvalidation over the experiments
#' Here each fold contains all the
#'  measured node values for a cue combination (there are 10 cue combinations in this data). 
R2 = crossvalidate(cnodata_prep, model = model,type = "experiment", nfolds=10)


#' we pick the worst performing fold to see what the code is doing

worseFold = which.max(R2$cvScores$cvScore)
whichI = R2$foldid == worseFold
timeIndex = 2
CNOlist.sub = cnodata_prep
CNOlist.cv = cnodata_prep

CNOlist.sub@signals[[timeIndex]][whichI,] = NA
CNOlist.cv@signals[[timeIndex]][!whichI,] = NA
#' data used for training:
plot(CNOlist.sub)

#' CV data:
plot(CNOlist.cv)

#' The vrossvalidation data corresponds to 1 experimental condition, which is missing from the training set.  

# fit the model on the fitting data
fit = gaBinaryT1(CNOlist = CNOlist.sub, 
				 model=model, timeIndex = timeIndex, verbose = F)

print(fit$bScore)
print(fit$bString)

#' comparing to the bitString of the optimal model
print(opt$bString)

#' optimal model's score on the training set:
optScore = computeScoreT1(CNOlist = CNOlist.sub,model = model,bString = opt$bString)
print(optScore)

#' show the model fit on the training set:
S1 = cutAndPlot(CNOlist.sub, model, bStrings=list(fit$bString))


#' model predictions on the cross-validation set
S1 = cutAndPlot(CNOlist.cv, model, bStrings=list(fit$bString))


#' compute score for the cross-validation
cvScore = computeScoreT1(CNOlist = CNOlist.cv,model = model,bString = fit$bString)
print(cvScore)



#' ## 3-fold crossvalidation over the measured nodes
#' Here each fold contains a subset of measured nodes. 
R3 = crossvalidate(cnodata_prep,model = model,type = "observable",nfolds=3)

#' we pick fold 3 to see what the code is doing

whichI = R3$foldid == 3
timeIndex = 2
CNOlist.sub = cnodata_prep
CNOlist.cv = cnodata_prep

CNOlist.sub@signals[[timeIndex]][,whichI] = NA
CNOlist.cv@signals[[timeIndex]][,!whichI] = NA


#' data used for training:
plot(CNOlist.sub)

#' CV data:
plot(CNOlist.cv)

#' We can see that the training data and the crossvalidation data overlaps regarding time point 0, but 
#' they are totally distinct at time point 1. 

#' fit the model on the fitting data
fit = gaBinaryT1(CNOlist = CNOlist.sub, 
				 model=model, timeIndex = timeIndex, verbose = F)

print(fit$bScore)
print(fit$bString)

#' comparing to the bitString of the optimal model
print(opt$bString)

#' optimal model's score on the training set:
optScore = computeScoreT1(CNOlist = CNOlist.sub,model = model,bString = opt$bString)
print(optScore)

#' show the model fit on the training set:
S1 = cutAndPlot(CNOlist.sub, model, bStrings=list(fit$bString))


#' model predictions on the cross-validation set
S1 = cutAndPlot(CNOlist.cv, model, bStrings=list(fit$bString))


#' compute score for the cross-validation
cvScore = computeScoreT1(CNOlist = CNOlist.cv,model = model,bString = fit$bString)

#' score for cross-validation:
print(cvScore)

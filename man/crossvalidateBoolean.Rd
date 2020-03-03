\name{crossvalidateBoolean}
\alias{crossvalidateBoolean}
\title{
k-fold crossvalidation for Boolean model.
}
\description{
Cross-validation analysis for the boolean case.
}
\usage{
crossvalidateBoolean = function(CNOlist,model,nfolds=10,foldid=NULL, 
                                type=c('datapoint','experiment','observable'),timeIndex = 2,parallel=FALSE, ...)
}
\arguments{
  \item{CNOlist}{
a CNOlist on which the score is based (based on valueSignals[[2]], i.e. data at time 1).
}
  \item{model}{
a model structure, as created by \code{readSIF}, normally pre-processed but that is not a requirement of this function.
}
  \item{nfolds}{
number of folds - default is 10. Although nfolds can be as large as the sample size (leave-one-out CV), it is not recommended for large datasets.
}
  \item{foldid}{
an optional vector of values between `1` and `nfold` identifying what fold each observation is in. If supplied, `nfold` can be missing.
}
  \item{type}{
define the way to do the crossvalidation. The default is type="datapoint"`, which assigns the data randomly into folds. The option `type="experiment"` uses whole experiments for crossvalidation (all data corresponding to a cue combination). The `type=observable` uses the subset of nodes across all experiments for crossvalidation.  
}
  \item{timeIndex}{
the index of the time point to optimize. Must be greater or equal to 2 (1 corresponds to time=0). Must be less than the number of time points. Default is 2.
}
  \item{parallel}{
verbose parameter, indicating wheter to parallelize the cross-validation procedure or not (default set to FALSE).
}


}
\details{
Does a k-fold cross-validation for Boolean CellNOpt models. In k-iterations a  fraction of the data is eliminated from the CNOlist. The model is trained on the  remaining data and then the model predicts the held-out data. Then the prediction accuracy is reported for each iteration.
}
\value{
This function returns a list with elements:
\item{cvScores}{cross-validation scores}
\item{fitScores}{fitting scores}
\item{bStrings}{the optimal bit-string list for each run}
\item{crossvalidate.call}{echo of the function which was called}
\item{foldid}{the fold id's}
}

\author{
A. Gabor, E. Gjerga
}
\examples{


data("ToyModel", package="CellNOptR")
data("CNOlistToy", package="CellNOptR")
pknmodel = ToyModel
cnodata = CNOlist(CNOlistToy)

# original and preprocessed network 
plotModel(pknmodel,cnodata)
model = preprocessing(data = cnodata,
					model = pknmodel,
					compression = TRUE,
					expansion = TRUE)
plotModel(model,cnodata)

# original CNOlist contains many timepoints, we use only a subset
plot(cnodata)
selectedTime = c(0,10)
cnodata_prep = cutCNOlist(cnodata,
                          model = model,
                          cutTimeIndices = which(!cnodata@timepoints \%in\% selectedTime))

plot(cnodata_prep)

# optimise and show results
opt = gaBinaryT1(CNOlist = cnodata_prep,model = model,verbose = FALSE)

# 10-fold crossvalidation using T1 data
# We use only T1 data for crossvalidation, because data in the T0 matrix is not independent.
# All rows of data in T0 describes the basal condition.

# Crossvalidation produce some text in the command window:  
\dontrun{
library(doParallel)
registerDoParallel(cores=3)
R=crossvalidateBoolean(CNOlist = cnodata_prep,
						model = model,
						type = "datapoint",
						nfolds = 10,
						parallel = TRUE)

}


}

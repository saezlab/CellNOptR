#' # adding permanentInhibitors and permanentStimuli functions
#' author: Attila Gabor  
#' date: 19.06.2018
#' 
#' The idea is to add permanent Cues, which can be turned on for time point 0  
#' simulations. 
#' For example, if we have an overexpressed experiment and we apply the stimuli, 
#' the overexpressed node should be activated already in T0. 
#' 
#' We modified the following files in CellNOptR
#' - CNOlist.R
#' - makeCNOlist.R
#' - computeScoreT1.R
#' - simulatorT0 (plotting)
load_all()
library(CellNOptR)

#' ## MacNamara case study:
pknmodel = readSIF("caseStudies/PKN-ToyPB.sif")
cnodata = CNOlist("caseStudies/MD-ToyPB.csv")
model = preprocessing(data = cnodata,model = pknmodel,compression = T,expansion = T)
plotModel(model,cnodata)
selectedTime = c(0,10)
cnodata_prep = cutCNOlist(cnodata,model = model,cutTimeIndices = which(!cnodata@timepoints %in% selectedTime))

#' optimise and show results
opt = gaBinaryT1(CNOlist = cnodata_prep,model = model,verbose = T)

plotModel(model,cnodata_prep,bString = opt$bString)
cutAndPlot(cnodata_prep,model,list(opt$bString))


#' ### Testing changes in computeScoreT1
bString = rep(1,length(model$reacID))
bString = opt$stringsTol[1,]

CNOlist = cnodata_prep
simList = prep4sim(model)
indexList = indexFinder(cnodata_prep, model)

bs = as.logical(bString)
modelCut <- list()
modelCut$interMat <- model$interMat[, bs]
modelCut$reacID <- model$reacID[bs]
modelCut$namesSpecies <- model$namesSpecies


simListCut <- cutSimList(simList, bString)


# Compute the simulated results
nStimuli = length(indexList$stimulated)
nInhibitors <- length(indexList$inhibited)
nCond <- dim(CNOlist@stimuli)[1]
nReacs <- length(modelCut$reacID)
nSpecies <- length(model$namesSpecies) # this is correct. No need to get modelCut$namesSpecies 
nMaxInputs <- dim(simListCut$finalCube)[2]


# simList matrices. C code must handle the matrix indices carefully.
# This is faster than transforming into a vector as in the previous code.
finalCube = as.integer(simListCut$finalCube-1)
ixNeg = as.integer(simListCut$ixNeg)
ignoreCube = as.integer(simListCut$ignoreCube)
maxIx = as.integer(simListCut$maxIx-1)

# index. convertion from R to C indices convention.
indexSignals <- as.integer(indexList$signals-1)
indexStimuli <- as.integer(indexList$stimulated-1)
indexInhibitors <- as.integer(indexList$inhibited-1)
nSignals <- length(indexSignals)


# cnolist
valueInhibitors <- as.integer(CNOlist@inhibitors)
valueStimuli <- as.integer(CNOlist@stimuli)

simResultsT1 = .Call("simulatorT1", nStimuli, nInhibitors,
				   nCond, nReacs, nSpecies, nSignals, nMaxInputs,
				   finalCube, ixNeg, ignoreCube, maxIx,
				   indexSignals, indexStimuli, indexInhibitors, valueInhibitors,
				   valueStimuli, as.integer(1))
simResultsT1

simResultsT0 = .Call("simulatorT1", nStimuli, nInhibitors,
					 nCond, nReacs, nSpecies, nSignals, nMaxInputs,
					 finalCube, ixNeg, ignoreCube, maxIx,
					 indexSignals, indexStimuli, indexInhibitors, 
					 valueInhibitors, valueStimuli, as.integer(0))
simResultsT0


# now is it possible to call simulatorT1 with mod=1 (to avoid zero out all the stimuli) and simulate the condition T0?
# we compare the original simulation (mod=as.integer(0)) vs zero-ing the stimulus matrix with mod=as.integer(1). 

valueInhibitorsT0 <- as.integer(matrix(0,nrow=nrow(CNOlist@inhibitors),ncol = ncol(CNOlist@inhibitors)))
valueStimuliT0 <-  as.integer(matrix(0,nrow=nrow(CNOlist@stimuli),ncol = ncol(CNOlist@stimuli)))

simResultsT0B = .Call("simulatorT1", nStimuli, nInhibitors,
					 nCond, nReacs, nSpecies, nSignals, nMaxInputs,
					 finalCube, ixNeg, ignoreCube, maxIx,
					 indexSignals, indexStimuli, indexInhibitors, 
					 valueInhibitorsT0, valueStimuliT0, as.integer(1))
simResultsT0B


# Now change EGF on in some experiments 
CNOlist@cues
CNOlist@stimuli
CNOlist@permanentStimuli[c(2,5,8),1]=1
CNOlist@permanentInhibitors[c(5),1]=0

valueInhibitorsT0 <-  CNOlist@inhibitors
valueInhibitorsT0[CNOlist@permanentInhibitors==0] = 0
valueInhibitorsT0 <- as.integer(valueInhibitorsT0)

valueStimuliT0 <-  CNOlist@stimuli
valueStimuliT0[CNOlist@permanentStimuli==0] = 0
valueStimuliT0 <- as.integer(valueStimuliT0)

simResultsT0C = .Call("simulatorT1", nStimuli, nInhibitors,
					  nCond, nReacs, nSpecies, nSignals, nMaxInputs,
					  finalCube, ixNeg, ignoreCube, maxIx,
					  indexSignals, indexStimuli, indexInhibitors, 
					  valueInhibitorsT0, valueStimuliT0, as.integer(1))
colnames(simResultsT0C) = model$namesSpecies
simResultsT0C
plotModel(model,cnodata,bString = bString)

#' ## TESTING
## all experiment except 2, 5 and 8 should be identical
diff = sum(abs(simResultsT0C[-c(2,5,8),] - simResultsT0[-c(2,5,8),]))
diff==0

## Experiment 2
CNOlist@cues[2,]
# T0 and T1 for exp 2 should be identical
diff = sum(abs(simResultsT0C[c(2),] - simResultsT1[c(2),]))
diff==0

## Experiment 5
CNOlist@cues[5,]
# EGF is active at T0 and but Pi3K is not inhibited. Identical to the exp2 @ T1
diff = sum(abs(simResultsT0C[c(5),] - simResultsT1[c(5),]))
diff==0


CNOlist@permanentInhibitors[c(5),1]=1

valueInhibitorsT0 <-  CNOlist@inhibitors
valueInhibitorsT0[CNOlist@permanentInhibitors==0] = 0
valueInhibitorsT0 <- as.integer(valueInhibitorsT0)


simResultsT0C = .Call("simulatorT1", nStimuli, nInhibitors,
					  nCond, nReacs, nSpecies, nSignals, nMaxInputs,
					  finalCube, ixNeg, ignoreCube, maxIx,
					  indexSignals, indexStimuli, indexInhibitors, 
					  valueInhibitorsT0, valueStimuliT0, as.integer(1))
colnames(simResultsT0C) = model$namesSpecies

## Experiment 5
CNOlist@cues[5,]
# EGF is active at T0 and Pi3K inhibited. Identical to the exp5 @ T1
diff = sum(abs(simResultsT0C[c(5),] - simResultsT1[c(5),]))
diff==0



# ---- end of testing computeScoreT1

#' ## Testing plot function
# cutAndPlot calls cutAndPlotResultsT1


# original:

cutAndPlotResultsT1(model = model,CNOlist = cnodata_prep,bString = opt$bString)
# permanent cues
CNOlist = cnodata_prep
CNOlist@permanentStimuli[c(2,5,8),1]=1
CNOlist@permanentInhibitors[c(5),1]=1

cutAndPlotResultsT1(model = model,CNOlist = CNOlist,bString = opt$bString)

plot(CNOlist)

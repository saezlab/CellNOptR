#' # Multi-model simulation for CellNOptR Boolean models
#' author: Attila Gabor  
#' date: 25.05.2018  

library(CellNOptR)
library(plyr)
library(dplyr)

#' ## MacNamara case study:
pknmodel = readSIF("caseStudies/PKN-ToyPB.sif")
cnodata = CNOlist("caseStudies/MD-ToyPB.csv")

#' preprocessed network 
model = preprocessing(data = cnodata,model = pknmodel,compression = T,expansion = T)
plotModel(model,cnodata)

#' original CNOlist contains many timepoints, we use only a subset
selectedTime = c(0,10)
cnodata_prep = cutCNOlist(cnodata,model = model,cutTimeIndices = which(!cnodata@timepoints %in% selectedTime))
plot(cnodata_prep)

#' optimise and show results
opt = gaBinaryT1(CNOlist = cnodata_prep,model = model,verbose = F)

#' optimal bitstring and score
print(opt$bString)
print(opt$bScore)

S = cutAndPlot(CNOlist = cnodata_prep,model = model,bStrings =list(opt$bString))

plotModel(model = model,CNOlist = cnodata_prep,bString = opt$bString)
plotModel(model = model,CNOlist = cnodata_prep,bString = colMeans(opt$stringsTol))



#' Distribution of scores in the tolerance
# compute score for all solutions in the Tolerance Set without the penalty

scores = apply(opt$stringsTol, MARGIN = 1, function(x){
	optScore = computeScoreT1(CNOlist = cnodata_prep,model = model,bString = x,sizeFac = 0,NAFac = 0)
} )

hist(scores)

#' from the histogram we can see at least 3 groups of models, which has different prediction score.

#' ## Visualisation of mixed model predictions  
#' keypoints:  
#' - do a plot similar to cutAndPlot, which shows the simulation/model fit of _multiple_ models/bitStrings  
#' - we dont want to show that CellNOpt predicts, let's say .7 activity of a node
#' - we want to show that there is an uncertainty in the prediction of certain nodes

# this will generate a matrix, where each column contains a simulation of a bitString for T0 and T1
ST1 = apply(opt$stringsTol, MARGIN = 1, function(x){
	s = simulateTN(CNOlist=cnodata_prep, model = model, bStrings = list(x))
	as.numeric(s)
} )

ST0 = apply(opt$stringsTol, MARGIN = 1, function(x){
	modelCut = cutModel(model, x)
	
	simList = prep4sim(model)
	simListCut <- cutSimList(simList, x)
	
	indexList=indexFinder(cnodata_prep,model = model)
	
	s = simulatorT0(CNOlist=cnodata_prep, model = modelCut, simList = simListCut,indexList = indexList)
	as.numeric(s)
} )

# mean and sd simulations
meanSimT0 = matrix(colMeans(t(ST0)),nrow = 10)
meanSimT1 = matrix(colMeans(t(ST1)),nrow = 10)
sdSimT0 = matrix(sqrt(diag(var(t(ST0)))),nrow = 10)
sdSimT1 = matrix(sqrt(diag(var(t(ST1)))),nrow = 10)

colnames(meanSimT0) = model$namesSpecies
colnames(meanSimT1) = model$namesSpecies
colnames(sdSimT0) = model$namesSpecies
colnames(sdSimT1) = model$namesSpecies
# colnames(sdSimT0) = paste0("SD_",model$namesSpecies)
# colnames(sdSimT1) = paste0("SD_",model$namesSpecies)


# combine the simulated and measured data frame for ggplot
indexList = indexFinder(CNOlist = cnodata_prep, model = model)
ST1.df = data.frame(meanSimT1[,indexList$signals])
ST0.df = data.frame(meanSimT0[,indexList$signals])
sdST1.df = data.frame(sdSimT1[,indexList$signals])
sdST0.df = data.frame(sdSimT0[,indexList$signals])



ST0.df$exps = 1:nrow(ST0.df)
ST0.df$T = "t0"
ST1.df$exps = 1:nrow(ST1.df)
ST1.df$T = "t1"
sdST0.df$exps = 1:nrow(sdST0.df)
sdST0.df$T = "t0"
sdST1.df$exps = 1:nrow(sdST1.df)
sdST1.df$T = "t1"


mP = rbind(ST0.df,ST1.df)
mP$source="sim"
sdmP = rbind(sdST0.df,sdST1.df)
sdmP$source="sim"

# mean and sd data
MT0.df = as.data.frame(cnodata_prep@signals[[1]])
MT1.df = as.data.frame(cnodata_prep@signals[[2]])
sdMT0.df = as.data.frame(cnodata_prep@variances[[1]])
sdMT1.df = as.data.frame(cnodata_prep@variances[[2]])
# colnames(sdMT0.df) = paste0("SD_",colnames(sdMT0.df))
# colnames(sdMT1.df) = paste0("SD_",colnames(sdMT1.df))

MT0.df$T = "t0"
MT1.df$T = "t1"
MT0.df$exps = 1:nrow(MT0.df)
MT1.df$exps = 1:nrow(MT1.df)
sdMT0.df$T = "t0"
sdMT1.df$T = "t1"
sdMT0.df$exps = 1:nrow(sdMT0.df)
sdMT1.df$exps = 1:nrow(sdMT1.df)


mp2 = rbind(MT0.df,MT1.df)
mp2$source="data"
sdmp2 = rbind(sdMT0.df,sdMT1.df)
sdmp2$source="data"


mP = rbind(mP,mp2)
sdmP = rbind(sdmP,sdmp2)



mPm = reshape2::melt(mP,id.vars=c("exps","T",'source'),variable.name="node")
sdmPm = reshape2::melt(sdmP,id.vars=c("exps","T",'source'),variable.name="node",value.name="SD")

mPm = merge.data.frame(mPm, sdmPm)

mPm2 = mPm
mPm2$value0 = 1-mPm2$value
mPm2$value1 = mPm2$value  # value=0 means, there is no model predicting 1. Tehrefore value1 = 0
mPm2$value=NULL
mPm2 = reshape2::melt(mPm2,id.vars=c('exps','T','node','source',"SD"))
#ggplot(mPm2) + geom_col(aes(x=T,y=value,fill=variable),col="black",position = 'dodge') + facet_grid(exps~node) + theme_bw() + scale_fill_manual(values=c('white','grey50'))


#' ### option 0: plot Mean and SD
#' 
ggplot(filter(mPm2,variable=="value1"), aes(x=T,y=value,col=source) ) +
	geom_point() +
	geom_line(aes(group=source)) +
	geom_errorbar(aes(ymin=value-SD, ymax=value+SD),width=0.2,size=1) + 
	#geom_crossbar(aes(ymin=value-SD, ymax=value+SD)) + 
	#geom_pointrange(aes(ymin=value-SD, ymax=value+SD),size=2) + 
	facet_grid(exps~node) + theme_bw() 

ggplot(filter(mPm2,variable=="value1"), aes(x=T,y=value,col=source) ) +
	geom_point() +
	geom_line(aes(group=source)) +
	#geom_errorbar(aes(ymin=value-SD, ymax=value+SD),width=0.2,size=1) + 
	geom_crossbar(aes(ymin=value-SD, ymax=value+SD),width=0.3) + 
	#geom_pointrange(aes(ymin=value-SD, ymax=value+SD),size=2) + 
	facet_grid(exps~node) + theme_bw() 

ggplot() + 
	geom_point(data=filter(mPm2,source=='sim',variable=="value1"),aes(x=T,y=value),col="darkgreen") +
	geom_errorbar(data=filter(mPm2,source=='sim',variable=="value1"),aes(x=T,ymin=value-SD, ymax=value+SD),col="darkgreen",width=.3) +
	geom_line(data=filter(mPm2,source=='sim',variable=="value1"),aes(x=T,y=value,group=node),col="darkgreen") +
	geom_point(data=filter(mPm2,source=='sim',variable=="value1"),aes(x=T,y=1,size=ifelse(value==0, NA, value)),col="darkgreen",stroke = 0, shape = 16) +
	geom_point(data=filter(mPm2,source=='sim',variable=="value1"),aes(x=T,y=0,size=ifelse((1-value)==0, NA, 1-value)),col="darkgreen",stroke = 0, shape = 16) +
	geom_point(data=filter(mPm2,source=='data',variable=="value1"),aes(x=T,y=value,size=1),col="tomato",stroke = 0, shape = 16) +
	geom_errorbar(data=filter(mPm2,source=='data',variable=="value1"),aes(x=T,ymin=value-SD, ymax=value+SD),col="tomato",width=.3) +
	geom_line(data=filter(mPm2,source=='data',variable=="value1"),aes(x=T,y=value,group=node),col="tomato") +
	facet_grid(exps~node) + theme_bw() + scale_size_area(max_size=5)


#' ### option 1: highlighting uncertainties
#' I mark the cases where the model predictions have some uncertainty.  
#' Problem:
#' - we cannot see the missfit in some nodes (ap1 in exp5/6 )
#' 
ggplot() + geom_col(data=filter(mPm2,source=='sim'),aes(x=T,y=value,fill=value,group=variable),col="black",position = 'dodge') +
	geom_point(data=filter(mPm2,source=='data',variable=="value1"),aes(x=T,y=value),col="black") +
	geom_line(data=filter(mPm2,source=='data',variable=="value1"),aes(x=T,y=value,group=node),col="black") +
	facet_grid(exps~node) + theme_bw() + scale_fill_gradient2(low="grey80",mid="red",high = 'grey80',midpoint = 0.5)


#' ### option 1b:  highlight missfit
#' Melanie's suggestion  
#' prediction uncertainty is shown by the barplots  
#' Needs update: two bars are colored the same way now.  
mPm3 = merge(filter(mPm2,source=='sim'),filter(mPm2,source=='data'),by = c('exps',"T", "node",'variable'))
mPm3$diff = abs(mPm3$value.x - mPm3$value.y)
ggplot() + geom_col(data=mPm3,aes(x=T,y=value.x,fill=diff,group=variable),col="black",position = 'dodge') +
	geom_point(data=filter(mPm3,variable=="value1"),aes(x=T,y=value.y),col="black") +
	geom_line(data=filter(mPm3,variable=="value1"),aes(x=T,y=value.y,group=node),col="black") +
	facet_grid(exps~node) + theme_bw() + scale_fill_gradient2(low="lightgreen",mid="white",high = 'firebrick',midpoint = 0.2)



#' ### option 2: highlight uncertainty using stacked bars
#' same issues as above, compact view
ggplot() + geom_col(data=filter(mPm2,source=='sim'),aes(x=T,y=value,fill=value,group=variable),col="black") +
	geom_point(data=filter(mPm2,source=='data',variable=="value1"),aes(x=T,y=value),col="black") +
	geom_line(data=filter(mPm2,source=='data',variable=="value1"),aes(x=T,y=value,group=node),col="black") +
	facet_grid(exps~node) + theme_bw() + scale_fill_gradient2(low="grey80",mid="red",high = 'grey80',midpoint = 0.5)


#' ### option 2b: highlight missift using stacked bars
#' Melanie's idea again
#' same issues as above, compact view  
#' _Needs update: two bars are colored the same way now._
mPm4 = mPm3
mPm4$variable =  as.character(mPm4$variable)
mPm4$variable[ mPm4$variable=='value0'] = "0"
mPm4$variable[ mPm4$variable=='value1'] = 1
mPm4$variable = as.numeric(mPm4$variable)
mPm4$diff = abs(mPm4$variable - mPm4$value.y)
ggplot() + geom_col(data=mPm3,aes(x=T,y=value.x,fill=diff>0.5),col="black") +
	geom_point(data=filter(mPm3,variable=="value1"),aes(x=T,y=value.y),col="black") +
	geom_line(data=filter(mPm3,variable=="value1"),aes(x=T,y=value.y,group=node),col="black") +
	facet_grid(exps~node) + theme_bw() + scale_fill_manual(values = c('green','red'))



#' ### option 3:
#' _outdataed._
ggplot() + geom_col(data=filter(mPm2,source=='sim'),aes(x=T,y=value,fill=variable,group=variable),col="black") +
	geom_point(data=filter(mPm2,source=='data',variable=="value1"),aes(x=T,y=value),col="black") +
	geom_line(data=filter(mPm2,source=='data',variable=="value1"),aes(x=T,y=value,group=node),col="black") +
	facet_grid(exps~node) + theme_bw() + scale_fill_manual(values = c('white','gray50'))

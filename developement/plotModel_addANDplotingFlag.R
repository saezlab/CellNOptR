#' # Flag to plot the eliminated AND gates when users draw the Model
#' author: Attila Gabor  
#' date: 25.05.2018  
load_all()
library(CellNOptR)

#' ## MacNamara case study:
pknmodel = readSIF("../caseStudies/PKN-ToyPB.sif")
cnodata = CNOlist("../caseStudies/MD-ToyPB.csv")

#' original and preprocessed network 

model = preprocessing(data = cnodata,model = pknmodel,compression = T,expansion = T)
plotModel(model,cnodata)

#' original CNOlist contains many timepoints, we use only a subset
selectedTime = c(0,10)
cnodata_prep = cutCNOlist(cnodata,model = model,cutTimeIndices = which(!cnodata@timepoints %in% selectedTime))

#' optimise and show results
opt = gaBinaryT1(CNOlist = cnodata_prep,model = model,verbose = F)

#' note that the edges corresponding to the AND gates are not presented
plotModel(model = model,CNOlist = cnodata_prep,bString = opt$bString)


#' 
debug(plotModel)
plotModel(model = model,CNOlist = cnodata_prep,bString = colMeans(opt$stringsTol),)

plotModel(model = model,CNOlist = cnodata_prep,bString = opt$bString,removeEmptyAnds = T)
plotModel(model = model,CNOlist = cnodata_prep,bString = opt$bString,removeEmptyAnds = F)

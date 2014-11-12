
# First, let us do the analysis manually, from scratch
# ----------------------------------------------------
#  
library(CellNOptR)
dir.create("CNOR_analysis")
setwd("CNOR_analysis")
cpfile<-dir(system.file("ToyModel",package="CellNOptR"),full=TRUE)
file.copy(from=cpfile,to=getwd(),overwrite=TRUE)
dataToy<-readMIDAS("ToyDataMMB.csv")
dataToy<-readMIDAS("ToyDataMMB.csv", verbose=FALSE)
CNOlistToy<-makeCNOlist(dataToy,subfield=FALSE)
data(CNOlistToy,package="CellNOptR")
plotCNOlist(CNOlistToy)
plotCNOlistPDF(CNOlist=CNOlistToy,filename="ToyModelGraph.pdf")
ToyModel<-readSIF("ToyPKNMMB.sif")
data(ToyModel,package="CellNOptR")
checkSignals(CNOlistToy,ToyModel)
ToyNCNOcutCompExp <- preprocessing(CNOlistToy, ToyModel, expansion=TRUE, compression=TRUE, cutNONC=TRUE, verbose=FALSE)
resECNOlistToy<-residualError(CNOlistToy)
initBstring<-rep(1,length(ToyNCNOcutCompExp$reacID))

ToyT1opt<-gaBinaryT1(CNOlist=CNOlistToy, model=ToyNCNOcutCompExp,initBstring=initBstring, verbose=FALSE)

cutAndPlot(model=ToyNCNOcutCompExp, bStrings=list(ToyT1opt$bString),CNOlist=CNOlistToy, plotPDF=TRUE)

plotFit(optRes=ToyT1opt)
pdf("evolFitToyT1.pdf")
plotFit(optRes=ToyT1opt)
dev.off()


writeScaffold(modelComprExpanded=ToyNCNOcutCompExp,optimResT1=ToyT1opt,optimResT2=NA, modelOriginal=ToyModel,CNOlist=CNOlistToy)

writeNetwork(modelOriginal=ToyModel,modelComprExpanded=ToyNCNOcutCompExp,optimResT1=ToyT1opt,optimResT2=NA,CNOlist=CNOlistToy)

namesFilesToy<-list(dataPlot="ToyModelGraph.pdf",evolFitT1="evolFitToyT1.pdf",evolFitT2=NA,simResultsT1="SimResultsT1_1.pdf",simResultsT2=NA,scaffold="Scaffold.sif",scaffoldDot="Scaffold.dot",tscaffold="TimesScaffold.EA",wscaffold="weightsScaffold.EA",PKN="PKN.sif",PKNdot="PKN.dot",wPKN="TimesPKN.EA",nPKN="nodesPKN.NA")

writeReport(modelOriginal=ToyModel,modelOpt=ToyNCNOcutCompExp,optimResT1=ToyT1opt,optimResT2=NA,CNOlist=CNOlistToy,directory="testToy",namesFiles=namesFilesToy,namesData=list(CNOlist="Toy",model="ToyModel"),resE=resECNOlistToy)

# The one step version --------------------------------
dataToy<-readMIDAS("ToyDataMMB.csv")
CNOlistToy<-makeCNOlist(dataToy,subfield=FALSE)
ToyModel<-readSIF("ToyPKNMMB.sif")

res <- CNORwrap(paramsList=NA, name="Toy", namesData=list(CNOlist="ToyData",model="ToyModel"),data=CNOlistToy, model=ToyModel)

pList = defaultParameters(CNOlistToy, ToyModel)
res <- CNORwrap( paramsList=pList,    name="Toy1Step",    namesData=list(CNOlist="ToyData",model="ToyModel"),    data=CNOlistToy,    model=ToyModel)

# real example -----------------------------
cpfile<-dir(system.file("DREAMModel",package="CellNOptR"),full=TRUE)
file.copy(from=cpfile,to=getwd(),overwrite=TRUE)
data(CNOlistDREAM,package="CellNOptR")
data(DreamModel,package="CellNOptR")
model = preprocessing(CNOlistDREAM, DreamModel, verbose=FALSE)
res = gaBinaryT1(CNOlistDREAM, model, verbose=FALSE, maxTime=10)
cutAndPlot(CNOlistDREAM, model, bStrings=list(res$bString), plotPDF=TRUE, plotParams=list(maxrow=25, margin=0.1, width=20, height=20))


# 2 time points -------------------------------
data(CNOlistToy2,package="CellNOptR")
data(ToyModel2,package="CellNOptR")
plotCNOlist(CNOlistToy2)
plotCNOlistPDF(CNOlist=CNOlistToy2,filename="ToyModelGraphT2.pdf")
ToyNCNOcutCompExp2 = preprocessing(CNOlistToy2,ToyModel2, verbose=FALSE)
initBstring2<-rep(1,length(ToyNCNOcutCompExp2$reacID))
ToyT1opt2<-gaBinaryT1(CNOlist=CNOlistToy2, model=ToyNCNOcutCompExp2, initBstring=initBstring2, stallGenMax=10, maxTime=60, verbose=FALSE)
cutAndPlot(model=ToyNCNOcutCompExp2, bStrings=list(ToyT1opt2$bString), CNOlist=CNOlistToy2, plotPDF=TRUE)


cutAndPlot(model=ToyNCNOcutCompExp2, bStrings=list(ToyT1opt2$bString), CNOlist=CNOlistToy2,  plotPDF=TRUE)
pdf("evolFitToy2T1.pdf")
plotFit(optRes=ToyT1opt2)
dev.off()
plotFit(optRes=ToyT1opt2)


ToyT1opt2T2<-gaBinaryTN(
    CNOlist=CNOlistToy2,
    model=ToyNCNOcutCompExp2,
    bStrings=list(ToyT1opt2$bString),
    stallGenMax=10,
    maxTime=60,
    verbose=FALSE)

cutAndPlot( model=ToyNCNOcutCompExp2, bStrings=list(ToyT1opt2$bString, ToyT1opt2T2$bString), CNOlist=CNOlistToy2,  plotPDF=TRUE)


pdf("evolFitToy2T2.pdf")
plotFit(optRes=ToyT1opt2T2)
dev.off()
plotFit(optRes=ToyT1opt2T2)

writeScaffold(
    modelComprExpanded=ToyNCNOcutCompExp2,
    optimResT1=ToyT1opt2,
    optimResT2=ToyT1opt2T2,
    modelOriginal=ToyModel2,
    CNOlist=CNOlistToy2)

writeNetwork(
    modelOriginal=ToyModel2,
    modelComprExpanded=ToyNCNOcutCompExp2,
    optimResT1=ToyT1opt2,
    optimResT2=ToyT1opt2T2,
    CNOlist=CNOlistToy2)

namesFilesToy<-list( dataPlot="ToyModelGraphT2.pdf", evolFitT1="evolFitToy2T1.pdf",   evolFitT2="evolFitToy2T2.pdf", simResultsT2="SimResultsTN.pdf",
    simResultsT1="SimResultsT1_1.pdf",    scaffold="Scaffold.sif",    scaffoldDot="Scaffold.dot",    tscaffold="TimesScaffold.EA",
    wscaffold="weightsScaffold.EA",    PKN="PKN.sif",    PKNdot="PKN.dot",   wPKN="TimesPKN.EA",    nPKN="nodesPKN.NA")

writeReport(modelOriginal=ToyModel2,    modelOpt=ToyNCNOcutCompExp2,optimResT1=ToyT1opt2,    optimResT2=ToyT1opt2T2,    CNOlist=CNOlistToy2,directory="testToy2",    namesFiles=namesFilesToy,namesData=list(CNOlist="ToyModified4T2",model="ToyModified4T2"))



library(CellNOptR)
dir.create("CNOR_analysis")
setwd("CNOR_analysis")
#If loading the data from the MIDAS file and the model from the sif file, type:
#cpfile<-dir(system.file("ToyModel",package="CellNOptR"),full=TRUE)
#file.copy(from=cpfile,to=getwd(),overwrite=TRUE)
#dataToy<-readMIDAS(MIDASfile='ToyDataMMB.csv')
#CNOlistToy<-makeCNOlist(dataset=dataToy,subfield=FALSE)
#ToyModel<-readSif(sifFile="ToyPKNMMB.sif")
#If loading data and model directly from the package:
data(CNOlistToy,package="CellNOptR")
data(ToyModel,package="CellNOptR")
CNOlistToy
plotCNOlist(CNOlistToy)
plotCNOlistPDF(CNOlist=CNOlistToy,filename="ToyModelGraph.pdf")
checkSignals(CNOlistToy,ToyModel)
indicesToy<-indexFinder(CNOlistToy,ToyModel,verbose=TRUE)
ToyNCNOindices<-findNONC(ToyModel,indicesToy,verbose=TRUE)
ToyNCNOcut<-cutNONC(ToyModel,ToyNCNOindices)
indicesToyNCNOcut<-indexFinder(CNOlistToy,ToyNCNOcut)
ToyNCNOcutComp<-compressModel(ToyNCNOcut,indicesToyNCNOcut)
indicesToyNCNOcutComp<-indexFinder(CNOlistToy,ToyNCNOcutComp)
ToyNCNOcutCompExp<-expandGates(ToyNCNOcutComp)
resECNOlistToy<-residualError(CNOlistToy)
ToyFields4Sim<-prep4Sim(ToyNCNOcutCompExp)
initBstring<-rep(1,length(ToyNCNOcutCompExp$reacID))
ToyT1opt<-gaBinaryT1(CNOlist=CNOlistToy,Model=ToyNCNOcutCompExp,SimList=ToyFields4Sim,indexList=indicesToyNCNOcutComp,initBstring=initBstring,verbose=TRUE)
cutAndPlotResultsT1(Model=ToyNCNOcutCompExp,bString=ToyT1opt$bString,SimList=ToyFields4Sim,CNOlist=CNOlistToy,indexList=indicesToyNCNOcutComp,plotPDF=TRUE)
plotFit(OptRes=ToyT1opt)
cutAndPlotResultsT1(Model=ToyNCNOcutCompExp,bString=ToyT1opt$bString,SimList=ToyFields4Sim,CNOlist=CNOlistToy,indexList=indicesToyNCNOcutComp,plotPDF=TRUE)
pdf("evolFitToyT1.pdf")
plotFit(OptRes=ToyT1opt)
dev.off()
writeScaffold(ModelComprExpanded=ToyNCNOcutCompExp,optimResT1=ToyT1opt,optimResT2=NA,ModelOriginal=ToyModel,CNOlist=CNOlistToy)
writeNetwork(ModelOriginal=ToyModel,ModelComprExpanded=ToyNCNOcutCompExp,optimResT1=ToyT1opt,optimResT2=NA,CNOlist=CNOlistToy)
namesFilesToy<-list(dataPlot="ToyModelGraph.pdf",evolFit1="evolFitToyT1.pdf",evolFit2=NA,SimResults1="ToyNCNOcutCompExpSimResultsT1.pdf",SimResults2=NA,Scaffold="ToyNCNOcutCompExpScaffold.sif",ScaffoldDot="ModelModelComprExpandedScaffold.dot",tscaffold="ToyNCNOcutCompExpTimesScaffold.EA",wscaffold="ToyNCNOcutCompExpweightsScaffold.EA",PKN="ToyModelPKN.sif",PKNdot="ToyModelPKN.dot",wPKN="ToyModelTimesPKN.EA",nPKN="ToyModelnodesPKN.NA")
writeReport(ModelOriginal=ToyModel,ModelOpt=ToyNCNOcutCompExp,optimResT1=ToyT1opt,optimResT2=NA,CNOlist=CNOlistToy,directory="testToy",namesFiles=namesFilesToy,namesData=list(CNOlist="Toy",Model="ToyModel"),resE=resECNOlistToy)
################################################
############The one step version################
library(CellNOptR)
dir.create("CNOR_analysis")
setwd("CNOR_analysis")
#If loading the data from the MIDAS file and the model from the sif file, type:
#cpfile<-dir(system.file("ToyModel",package="CellNOptR"),full=TRUE)
#file.copy(from=cpfile,to=getwd(),overwrite=TRUE)
#dataToy<-readMIDAS(MIDASfile='ToyDataMMB.csv')
#CNOlistToy<-makeCNOlist(dataset=dataToy,subfield=FALSE)
#ToyModel<-readSif(sifFile="ToyPKNMMB.sif")
#If loading data and model directly from the package:
data(CNOlistToy,package="CellNOptR")
data(ToyModel,package="CellNOptR")
CNORwrap(paramsList=NA,Name="Toy",NamesData=list(CNOlist="ToyData",Model="ToyModel"),Data=CNOlistToy,Model=ToyModel)
#version 2
pList<-list(Data=CNOlistToy,Model=ToyModel,sizeFac = 1e-04, NAFac = 1, PopSize = 50, Pmutation = 0.5, MaxTime = 60, maxGens = 500, StallGenMax = 100, SelPress = 1.2, elitism = 5, RelTol = 0.1,verbose=TRUE)
CNORwrap(paramsList=pList,Name="Toy",NamesData=list(CNOlist="ToyData",Model="ToyModel"),Data=NA,Model=NA)
######################################################
############The DREAM data and network################
library(CellNOptR)
dir.create("CNOR_analysis")
setwd("CNOR_analysis")
#If loading the data from the MIDAS file and the model from the sif file, type:
#cpfile<-dir(system.file("DREAMModel",package="CellNOptR"),full=TRUE)
#file.copy(from=cpfile,to=getwd(),overwrite=TRUE)
#dreamData<-readMIDAS(MIDASfile="LiverDataDREAMMIDAS.csv")
#CNOlistDREAM<-makeCNOlist(dataset=dreamData,subfield=FALSE)
#DreamModel<-readSif(sifFile="LiverPKNDREAM.sif")
#If loading data and model directly from the package:
data(CNOlistDREAM,package="CellNOptR")
data(DreamModel,package="CellNOptR")
plotCNOlist(CNOlistDREAM)
plotCNOlistPDF(CNOlist=CNOlistDREAM,filename="DREAMdataGraph.pdf")
checkSignals(CNOlistDREAM,DreamModel)
indicesDream<-indexFinder(CNOlistDREAM,DreamModel,verbose=TRUE)
#Find the index of the non-obs/non-ctrl
DreamNCNOindices<-findNONC(DreamModel,indicesDream,verbose=TRUE)
#cut the nonc off the model
DreamNCNOcut<-cutNONC(DreamModel,DreamNCNOindices)
#recompute the indices of the species
indicesDreamNCNOcut<-indexFinder(CNOlistDREAM,DreamNCNOcut)
#compress the model
DreamNCNOcutComp<-compressModel(DreamNCNOcut,indicesDreamNCNOcut)
#Recompute the indices of the species
indicesDreamNCNOcutComp<-indexFinder(CNOlistDREAM,DreamNCNOcutComp)
#Expand the gates
DreamNCNOcutCompExp<-expandGates(DreamNCNOcutComp)
resECNOdream<-residualError(CNOlistDREAM)
#Prepare for simulation
DreamFields4Sim<-prep4Sim(DreamNCNOcutCompExp)
#Optimisation
initBstring<-rep(1,length(DreamNCNOcutCompExp$reacID))
DreamT1opt<-gaBinaryT1(CNOlist=CNOlistDREAM,Model=DreamNCNOcutCompExp,SimList=DreamFields4Sim,indexList=indicesDreamNCNOcutComp,initBstring=initBstring,verbose=TRUE,MaxTime=600)
#Plot simulated and experimental results
cutAndPlotResultsT1(Model=DreamNCNOcutCompExp,bString=DreamT1opt$bString,SimList=DreamFields4Sim,CNOlist=CNOlistDREAM,indexList=indicesDreamNCNOcutComp,plotPDF=TRUE)
#Plot the evolution of fit
pdf("evolFitDreamT1.pdf")
plotFit(OptRes=DreamT1opt)
dev.off()
plotFit(OptRes=DreamT1opt)
writeScaffold(ModelComprExpanded=DreamNCNOcutCompExp,optimResT1=DreamT1opt,optimResT2=NA,ModelOriginal=DreamModel,CNOlist=CNOlistDREAM)
writeNetwork(ModelOriginal=DreamModel,ModelComprExpanded=DreamNCNOcutCompExp,optimResT1=DreamT1opt,optimResT2=NA,CNOlist=CNOlistDREAM)
#Write the report
namesFilesDream<-list(dataPlot="DREAMdataGraph.pdf",evolFit1="evolFitDreamT1.pdf",evolFit2=NA,SimResults1="DreamNCNOcutCompExpSimResultsT1.pdf",SimResults2=NA,Scaffold="DreamNCNOcutCompExpScaffold.sif",ScaffoldDot="ModelModelComprExpandedScaffold.dot",tscaffold="DreamNCNOcutCompExpTimesScaffold.EA",wscaffold="DreamNCNOcutCompExpweightsScaffold.EA",PKN="DreamModelPKN.sif",PKNdot="DreamModelPKN.dot",wPKN="DreamModelTimesPKN.EA",nPKN="DreamModelnodesPKN.NA")
writeReport(ModelOriginal=DreamModel,ModelOpt=DreamNCNOcutCompExp,optimResT1=DreamT1opt,optimResT2=NA,CNOlist=CNOlistDREAM,directory="testDREAM",namesFiles=namesFilesDream,namesData=list(CNOlist="Dream",Model="DreamModel"),resE=resECNOdream)
################################################
############The 2 time points################
library(CellNOptR)
dir.create("CNOR_analysis")
setwd("CNOR_analysis")
######
cpfile<-dir(system.file("ToyModel",package="CellNOptR"),full=TRUE)
file.copy(from=cpfile,to=getwd(),overwrite=TRUE)
dataToy<-readMIDAS(MIDASfile='ToyDataMMB.csv')
CNOlistToy<-makeCNOlist(dataset=dataToy,subfield=FALSE)
CNOlistToy
#Transform data for multiple time points
CNOlistToy2<-CNOlistToy
CNOlistToy2$valueSignals[[3]]<-CNOlistToy2$valueSignals[[2]]
CNOlistToy2$valueSignals[[3]][,6:7]<-0
CNOlistToy2$valueSignals[[2]][which(CNOlistToy2$valueSignals[[2]][,6] > 0),6]<-0.5
CNOlistToy2$valueSignals[[2]][which(CNOlistToy2$valueSignals[[2]][,7] > 0),7]<-0.77118
CNOlistToy2$timeSignals<-c(CNOlistToy2$timeSignals, 100)
#In this model I added a negative fedback between cJun and Jnk (!cJun=Jnk)
#this is the model to use with the data CNOlistToy2
ToyModel2<-readSif(sifFile="ToyModelMMB2.sif")
#####
data(CNOlistToy2,package="CellNOptR")
data(ToyModel2,package="CellNOptR")
#
plotCNOlist(CNOlistToy2)
plotCNOlistPDF(CNOlist=CNOlistToy2,filename="ToyModelGraphT2.pdf")
checkSignals(CNOlistToy2,ToyModel2)
indexesToy2<-indexFinder(CNOlistToy2,ToyModel2,verbose=TRUE)
ToyNCNOindices2<-findNONC(ToyModel2,indexesToy2,verbose=TRUE)
ToyNCNOcut2<-cutNONC(ToyModel2,ToyNCNOindices2)
indexesToyNCNOcut2<-indexFinder(CNOlistToy2,ToyNCNOcut2)
ToyNCNOcutComp2<-compressModel(ToyNCNOcut2,indicesToyNCNOcut2)
indicesToyNCNOcutComp2<-indexFinder(CNOlistToy2,ToyNCNOcutComp2)
ToyNCNOcutCompExp2<-expandGates(ToyNCNOcutComp2)
resECNOlistToy2<-residualError(CNOlistToy2)
ToyFields4Sim2<-prep4Sim(ToyNCNOcutCompExp2)
initBstring2<-rep(1,length(ToyNCNOcutCompExp2$reacID))
ToyT1opt2<-gaBinaryT1(CNOlist=CNOlistToy2,Model=ToyNCNOcutCompExp2,SimList=ToyFields4Sim2,indexList=indicesToyNCNOcutComp2,initBstring=initBstring2,MaxTime=180)
cutAndPlotResultsT1(Model=ToyNCNOcutCompExp2,bString=ToyT1opt2$bString,SimList=ToyFields4Sim2,CNOlist=CNOlistToy2,indexList=indicesToyNCNOcutComp2,plotPDF=TRUE)
pdf("evolFitToy2T1.pdf")
plotFit(OptRes=ToyT1opt2)
dev.off()
plotFit(OptRes=ToyT1opt2)
#Optimise T2
SimToyT12<-simulateT1(CNOlist=CNOlistToy2,Model=ToyNCNOcutCompExp2,bStringT1=ToyT1opt2$bString,SimList=ToyFields4Sim2,indexList=indicesToyNCNOcutComp2)
ToyT1opt2T2<-gaBinaryT2(CNOlist=CNOlistToy2,Model=ToyNCNOcutCompExp2,SimList=ToyFields4Sim2,indexList=indicesToyNCNOcutComp2,bStringT1=ToyT1opt2$bString,SimResT1=SimToyT12,MaxTime=180)
cutAndPlotResultsT2(Model=ToyNCNOcutCompExp2,bStringT1=ToyT1opt2$bString,bStringT2=ToyT1opt2T2$bString,SimList=ToyFields4Sim2,CNOlist=CNOlistToy2,indexList=indicesToyNCNOcutComp2,plotPDF=TRUE)
pdf("evolFitToy2T2.pdf")
plotFit(OptRes=ToyT1opt2T2)
dev.off()
plotFit(OptRes=ToyT1opt2T2)
writeScaffold(ModelComprExpanded=ToyNCNOcutCompExp2,optimResT1=ToyT1opt2,optimResT2=ToyT1opt2T2,ModelOriginal=ToyModel2,CNOlist=CNOlistToy2)
writeNetwork(ModelOriginal=ToyModel2,ModelComprExpanded=ToyNCNOcutCompExp2,optimResT1=ToyT1opt2,optimResT2=ToyT1opt2T2,CNOlist=CNOlistToy2)
namesFilesToy<-list(dataPlot="ToyModelGraphT2.pdf",evolFit1="evolFitToy2T1.pdf",evolFit2="evolFitToy2T2.pdf",SimResults2="ToyNCNOcutCompExp2SimResultsT1T2.pdf",SimResults1="ToyNCNOcutCompExp2SimResultsT1.pdf",Scaffold="Scaffold.sif",ScaffoldDot="Scaffold.dot",tscaffold="TimesScaffold.EA",wscaffold="weightsScaffold.EA",PKN="PKN.sif",PKNdot="PKN.dot",wPKN="TimesPKN.EA",nPKN="nodesPKN.NA")
writeReport(ModelOriginal=ToyModel2,ModelOpt=ToyNCNOcutCompExp2,optimResT1=ToyT1opt2,optimResT2=ToyT1opt2T2,CNOlist=CNOlistToy2,directory="testToy2",namesFiles=namesFilesToy,namesData=list(CNOlist="ToyModified4T2",Model="ToyModified4T2"),resE=resECNOlistToy2)
################################################
############The 2 time points -- one step################
library(CellNOptR)
dir.create("CNOR_analysis")
setwd("CNOR_analysis")
data(CNOlistToy2,package="CellNOptR")
data(ToyModel2,package="CellNOptR")
pList<-list(Data=CNOlistToy2,Model=ToyModel2,sizeFac = 1e-04, NAFac = 1, PopSize = 50, Pmutation = 0.5, MaxTime = 60, maxGens = 500, StallGenMax = 100, SelPress = 1.2, elitism = 5, RelTol = 0.1,verbose=TRUE)
CNORwrap(paramsList=pList,Name="Toy",NamesData=list(CNOlist="ToyData2",Model="ToyModel2"),Data=NA,Model=NA,Time=2)


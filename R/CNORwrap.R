#This function is a wrapper around the whole CNOR analysis, it performs the following steps:
#1.Plot the CNOlist
#2.Checks data to model compatibility
#3.Find the indices, in the model, of the species that are inh/stim/sign
#4.Find the indices of the non-osb/non-contr
#5.Cut the nonc off the model
#6.Recompute the indices
#7.Compress the model
#8.Recompute the indices
#9.Expand the gates
#10.Compute the residual error
#11.Prepare for simulation
#12.Optimisation t1
#13.Plot simulated and experimental results
#14.Plot the evolution of fit
#15.Optimise t2 (not implemented in this version)
#16.Write the scaffold and PKN
#17.Write the report

CNORwrap<-function(paramsList,Data,Model,Name,NamesData,Time=1){

#if the paramsList is set to NA, it means that only default parameters have been provided
#so we are going to build the parameters list with those

	if(is.na(paramsList[1])){
		paramsList<-list()
		paramsList$Data<-Data
		paramsList$Model<-Model
		paramsList$sizeFac<-1e-04
		paramsList$NAFac<-1
		paramsList$PopSize<-50
		paramsList$Pmutation<-0.5
		paramsList$MaxTime<-60
		paramsList$maxGens<-500
		paramsList$StallGenMax<-100
		paramsList$SelPress<-1.2
		paramsList$elitism<-5
		paramsList$RelTol<-0.1
		paramsList$verbose<-FALSE
		}
		
#1.Plot the CNOlist
	plotCNOlist(paramsList$Data)
	plotCNOlistPDF(
		CNOlist=paramsList$Data,
		fileName=paste(Name,"DataPlot.pdf",sep="")
		)
		
#2.Checks data to model compatibility
	checkSignals(CNOlist=paramsList$Data,Model=paramsList$Model)
	
#3.Find the indices, in the model, of the species that are inh/stim/sign
	Indices<-indexFinder(
		CNOlist=paramsList$Data,
		Model=paramsList$Model,
		verbose=paramsList$verbose)
	
#4. Find the indices of the non-osb/non-contr	
	NCNOindices<-findNONC(
		Model=paramsList$Model,
		indexes=Indices,
		verbose=paramsList$verbose)
		
#5.Cut the nonc off the model
	NCNOcut<-cutNONC(Model=paramsList$Model, NONCindexes=NCNOindices)
	
#6.Recompute the indices
	IndicesNCNOcut<-indexFinder(CNOlist=paramsList$Data,Model=NCNOcut)
	
#7.Compress the model
	NCNOcutComp<-compressModel(Model=NCNOcut,indexes=IndicesNCNOcut)
	
#8.Recompute the indices
	IndicesNCNOcutComp<-indexFinder(CNOlist=paramsList$Data,Model=NCNOcutComp)
	
#9.Expand the gates	
	NCNOcutCompExp<-expandGates(Model=NCNOcutComp)
	
#10.Compute the residual error
	resE<-residualError(CNOlist=paramsList$Data)
	
#11.Prepare for simulation
	fields4Sim<-prep4Sim(Model=NCNOcutCompExp)
	
#12.Optimisation t1	
	initBstring<-rep(1,length(NCNOcutCompExp$reacID))
	T1opt<-gaBinaryT1(CNOlist=paramsList$Data,
		Model=NCNOcutCompExp,
		SimList=fields4Sim,
		indexList=IndicesNCNOcutComp,
		initBstring=initBstring,
		sizeFac=paramsList$sizeFac,
		NAFac=paramsList$NAFac,
		PopSize=paramsList$PopSize,
		Pmutation=paramsList$Pmutation,
		MaxTime=paramsList$MaxTime,
		maxGens=paramsList$maxGens,
		StallGenMax=paramsList$StallGenMax,
		SelPress=paramsList$SelPress,
		elitism=paramsList$elitism,
		RelTol=paramsList$RelTol,
		verbose=paramsList$verbose)
		
#13.Plot simulated and experimental results
	cutAndPlotResultsT1(
		Model=NCNOcutCompExp,
		bString=T1opt$bString,
		SimList=fields4Sim,
		CNOlist=paramsList$Data,
		indexList=IndicesNCNOcutComp,
		plotPDF=TRUE)
		
#14.Plot the evolution of fit
	pdf(paste(Name,"evolFitT1.pdf",sep=""))
	plotFit(OptRes=T1opt)
	dev.off()
	plotFit(OptRes=T1opt)
	
#15.Optimise t2	
#Not implemented in this version
	T2opt<-NA
	
#16.Write the scaffold and PKN
#and
#17.Write the report
	if(Time==2){
#not implemented in this version

		writeScaffold(
			ModelComprExpanded=NCNOcutCompExp,
			optimResT1=T1opt,
			optimResT2=T2opt,
			ModelOriginal=paramsList$Model,
			CNOlist=paramsList$Data)
		writeNetwork(
			ModelOriginal=paramsList$Model,
			ModelComprExpanded=NCNOcutCompExp,
			optimResT1=T1opt,
			optimResT2=T2opt,
			CNOlist=paramsList$Data)
		namesfiles<-list(
			dataPlot=paste(Name,"DataPlot.pdf",sep=""),
			evolFit1=paste(Name,"evolFitT1.pdf",sep=""),
			evolFit2=paste(Name,"evolFitT2.pdf",sep=""),
			SimResults2="NCNOcutCompExpSimResultsT1T2.pdf",
			SimResults1="NCNOcutCompExpSimResultsT1.pdf",
			Scaffold="Scaffold.sif",
			ScaffoldDot="Scaffold.dot",
			tscaffold="TimesScaffold.EA",
			wscaffold="weightsScaffold.EA",
			PKN="PKN.sif",
			PKNdot="PKN.dot",
			wPKN="TimesPKN.EA",
			nPKN="nodesPKN.NA")
		writeReport(
			ModelOriginal=paramsList$Model,
			ModelOpt=NCNOcutCompExp,
			optimResT1=T1opt,
			optimResT2=T2opt,
			CNOlist=paramsList$Data,
			directory=Name,
			namesFiles=namesfiles,
			namesData=NamesData,
			resE=resE)	
	
		}else{
		
			writeScaffold(
				ModelComprExpanded=NCNOcutCompExp,
				optimResT1=T1opt,
				optimResT2=NA,
				ModelOriginal=paramsList$Model,
				CNOlist=paramsList$Data)
			writeNetwork(
				ModelOriginal=paramsList$Model,
				ModelComprExpanded=NCNOcutCompExp,
				optimResT1=T1opt,
				optimResT2=NA,
				CNOlist=paramsList$Data)
			namesfiles<-list(
				dataPlot=paste(Name,"DataPlot.pdf",sep=""),
				evolFit1=paste(Name,"evolFitT1.pdf",sep=""),
				evolFit2=NA,SimResults2=NA,
				SimResults1="NCNOcutCompExpSimResultsT1.pdf",
				Scaffold="Scaffold.sif",
				ScaffoldDot="Scaffold.dot",
				tscaffold="TimesScaffold.EA",
				wscaffold="weightsScaffold.EA",
				PKN="PKN.sif",
				PKNdot="PKN.dot",
				wPKN="TimesPKN.EA",
				nPKN="nodesPKN.NA")
			writeReport(
				ModelOriginal=paramsList$Model,
				ModelOpt=NCNOcutCompExp,
				optimResT1=T1opt,
				optimResT2=NA,
				CNOlist=paramsList$Data,
				directory=Name,
				namesFiles=namesfiles,
				namesData=NamesData,
				resE=resE)
	
			}
}
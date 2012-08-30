#
#  This file is part of the CNO software
#
#  Copyright (c) 2011-2012 - EBI
#
#  File author(s): CNO developers (cno-dev@ebi.ac.uk)
#
#  Distributed under the GPLv2 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-2.0.html
#
#  CNO website: http://www.ebi.ac.uk/saezrodriguez/software.html
#
##############################################################################
# $Id: CNORwrap.R 2256 2012-08-29 15:49:30Z cokelaer $
#This function is a wrapper around the whole CNOR analysis, it performs the following steps:
#1.Plot the CNOlist
#2.Checks data to model compatibility
#3.Cut the nonc off the model
#4.Compress the model
#5.Expand the gates
#6.Compute the residual error
#7.Prepare for simulation
#8.Optimisation t1
#9.Plot simulated and experimental results
#10.Plot the evolution of fit
#11.Optimise t2 (not implemented in this version)
#12.Write the scaffold and PKN
#13.Write the report

CNORwrap<-function(paramsList=NA, data=NA, model=NA, name, namesData=NA, time=1,
compression=TRUE, expansion=TRUE, cutNONC=TRUE)
{

    # aliases
    Name = name
    NamesData = namesData
    Time = time

    # if paramsList empty, we will fill it with
    # default values and the Data and model provided
    if(is.na(paramsList[1])==TRUE){
        #print("paramList not provided")
        # Data must be provided
        # is.na raise warning, so let us use is.list
        if (is.na(data[1])==TRUE){
            stop("if paramsList not provided, Data must be provided")
        }
        # model must be provided
        if (is.na(model[1])==TRUE){
            stop("if paramsList not provided, model must be provided")
        }
        paramsList = defaultParameters(data, model)
    }

    # if paramsList is provided, and model and Data are NA, then we must find
    # Data and model in paramsList
    else if (is.na(paramsList[1])==FALSE){
        #print("paramsList  provided")

        if (is.na(data[1])==TRUE && is.na(paramsList$data[1])==TRUE){
            stop("Data must be provided either in paramsList or Data argument")
        }
        if (is.na(model[1])==TRUE && is.na(paramsList$model[1])==TRUE){
            stop("model must be provided either in paramsList or model argument")
        }

        # Data must be provided
        if (is.na(data[1])==FALSE){
            #print("Data provided")
            if (is.na(paramsList$data[1])==FALSE){
                warning("overwritting paramsList$data with provided Data")
            }
            paramsList$data = data
        }
        # model must be provided
        if (is.na(model[1])==FALSE){
            #print("model provided")
            if (is.na(paramsList$model[1])==FALSE){
                warning("overwritting paramsList$Model with provided Model")
            }
            paramsList$model = model
        }

    }

    if (is.list(NamesData)==FALSE){
        if (is.na(NamesData)==TRUE){
            NamesData <- list(
                CNOlist=paste(Name, "Data",sep=""),
               model=paste(Name,"Model", sep=""))
        }
        else{
            stop("NamesData must be a list or kept to the default values (NA)")
        }
    }


    #1.Plot the CNOlist
    plotCNOlist(paramsList$data)
    plotCNOlistPDF(
        CNOlist=paramsList$data,
        filename=paste(Name,"DataPlot.pdf",sep="")
        )

    #2. Checks data to model compatibility
    checkSignals(CNOlist=paramsList$data,model=paramsList$model)

    #3.Cut the nonc off the model
    #4.Compress the model
    #5.Expand the gates
    newmodel = preprocessing(paramsList$data, paramsList$model,
        compression=compression, expansion=expansion, cutNONC=cutNONC)
    NCNOcutCompExp <- newmodel

    #6.Compute the residual error

    #7.Prepare for simulation

    #8.Optimisation t1
    initBstring<-rep(1,length(NCNOcutCompExp$reacID))
    T1opt<-gaBinaryT1(CNOlist=paramsList$data,
        model=NCNOcutCompExp,
        initBstring=initBstring,
        sizeFac=paramsList$sizeFac,
        NAFac=paramsList$NAFac,
        popSize=paramsList$popSize,
        pMutation=paramsList$pMutation,
        maxTime=paramsList$maxTime,
        maxGens=paramsList$maxGens,
        stallGenMax=paramsList$stallGenMax,
        selPress=paramsList$selPress,
        elitism=paramsList$elitism,
        relTol=paramsList$relTol,
        verbose=paramsList$verbose)

    #9.Plot simulated and experimental results
    cutAndPlot(
        model=NCNOcutCompExp,
        bStrings=list(T1opt$bString),
        CNOlist=paramsList$data,
        plotPDF=TRUE)


    #10.Plot the evolution of fit
    pdf(paste(Name,"evolFitT1.pdf",sep=""))
    plotFit(optRes=T1opt)
    dev.off()
    plotFit(optRes=T1opt)

    #11.Optimise t2
    if(Time==2){

        T2opt<-gaBinaryT2(
            CNOlist=paramsList$data,
            model=NCNOcutCompExp,
            bStringT1=T1opt$bString,
            sizeFac=paramsList$sizeFac,
            NAFac=paramsList$NAFac,
            popSize=paramsList$popSize,
            pMutation=paramsList$pMutation,
            maxTime=paramsList$maxTime,
            maxGens=paramsList$maxGens,
            stallGenMax=paramsList$stallGenMax,
            selPress=paramsList$selPress,
            elitism=paramsList$elitism,
            relTol=paramsList$relTol,
            verbose=paramsList$verbose)

        cutAndPlot(model=NCNOcutCompExp,CNOlist=paramsList$data,
            bStrings=list(T1opt$bString, bStringT2=T2opt$bString), plotPDF=TRUE)

        pdf(paste(Name,"evolFitT2.pdf",sep=""))
        plotFit(optRes=T2opt)
        dev.off()
        plotFit(optRes=T2opt)

        }
    else{
        T2opt<-NA
    }
    #13.Write the scaffold and PKN
    #and
    #14.Write the report
    writeScaffold(
        modelComprExpanded=NCNOcutCompExp,
        optimResT1=T1opt,
        optimResT2=T2opt,
        modelOriginal=paramsList$model,
        CNOlist=paramsList$data)
    writeNetwork(
        modelOriginal=paramsList$model,
        modelComprExpanded=NCNOcutCompExp,
        optimResT1=T1opt,
        optimResT2=T2opt,
        CNOlist=paramsList$data)

    if(Time==2){

        namesfiles<-list(
            dataPlot=paste(Name,"DataPlot.pdf",sep=""),
            evolFitT1=paste(Name,"evolFitT1.pdf",sep=""),
            evolFitT2=paste(Name,"evolFitT2.pdf",sep=""),
            simResultsT2="SimResultsTN.pdf",
            simResultsT1="SimResultsT1_1.pdf",
            scaffold="Scaffold.sif",
            scaffoldDot="Scaffold.dot",
            tscaffold="TimesScaffold.EA",
            wscaffold="weightsScaffold.EA",
            PKN="PKN.sif",
            PKNdot="PKN.dot",
            wPKN="TimesPKN.EA",
            nPKN="nodesPKN.NA")

    }
    else{
        namesfiles<-list(
            dataPlot=paste(Name,"DataPlot.pdf",sep=""),
            evolFitT1=paste(Name,"evolFitT1.pdf",sep=""),
            evolFitT2=NA,
            simResultsT2=NA,
            simResultsT1="SimResultsT1_1.pdf",
            scaffold="Scaffold.sif",
            scaffoldDot="Scaffold.dot",
            tscaffold="TimesScaffold.EA",
            wscaffold="weightsScaffold.EA",
            PKN="PKN.sif",
            PKNdot="PKN.dot",
            wPKN="TimesPKN.EA",
            nPKN="nodesPKN.NA")
    }

    writeReport(
        modelOriginal=paramsList$model,
        modelOpt=NCNOcutCompExp,
        optimResT1=T1opt,
        optimResT2=T2opt,
        CNOlist=paramsList$data,
        directory=Name,
        namesFiles=namesfiles,
        namesData=NamesData)

}

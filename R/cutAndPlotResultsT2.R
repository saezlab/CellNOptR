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
# $Id$

cutAndPlotResultsT2 <-function(model, bStringT1, bStringT2, CNOlist, simList=NULL,
    indexList=NULL, plotPDF=FALSE, tag=NULL, 
    tPt=CNOlist$timeSignals[2:3], plotParams=list(maxrow=10))
{
    warning("cutAndPlotResultsT2 is a deprecated function. Use cutAndPlotResultsTN instead. ")

    # ideally, CNOlist must be first argument but this function is deprecated
    # anyway, so keep as it is for back compatibility.
    if (is.null(simList)==TRUE){
        simList = prep4sim(model)
    }
    if (is.null(indexList)==TRUE){
        indexList = indexFinder(CNOlist, model)
    }

    if ("maxrow" %in% names(plotParams) == FALSE){
        plotParams$maxrow = 10
    }   



    modelCut <- cutModel(model, bStringT1)
    simListCut <- cutSimList(simList, bStringT1)

    # t0
    Sim0 <- simulatorT0(CNOlist=CNOlist,model=modelCut,simList=simListCut,indexList=indexList)
    simResT0 <- as.matrix(Sim0[,indexList$signals])

    # t1
    # same as simulateTN followed by selection of indexList$signals
    # SimT1 = simulatorT1(CNOlist=CNOlist, model=modelCut, simList=simListCut, indexList=indexList)
    SimT1 = simulateTN(CNOlist, model, bStrings=list(bStringT1))
    simResT1 <- as.matrix(SimT1[,indexList$signals])

    SimT2 = simulateTN(CNOlist, model, bStrings=list(bStringT1, bStringT2))
    simResT2 <- as.matrix(SimT2[,indexList$signals])

    # put it all together
    simResults <- list(
        t0=simResT0,
        t1=simResT1,
        t2=simResT2)


    # if there is a lot of data, split up cnolist
    # make the max dimensions 10 x 10

    dim1 = dim(CNOlist$valueSignals[[1]])[1]
    dim2 = dim(CNOlist$valueSignals[[1]])[2]

    CNOlistSet = list()
    simResultsSet = list()

    if(dim1 > plotParams$maxrow) { #|| dim2 > 10) {

        par1 = ceiling(dim1/plotParams$maxrow)
        div1 = ceiling(dim1/par1)
        par2 = ceiling(dim2/plotParams$maxrow)
        div2 = ceiling(dim2/par2)

        count1 = 1
        for(a in 1:par1) {
            CNOdiv = CNOlist
            simDiv = simResults
            finalN = div1 * a
            if(finalN > dim1) {finalN = dim1}
            CNOdiv$valueCues = CNOdiv$valueCues[count1:finalN,]
            CNOdiv$valueStimuli = CNOdiv$valueStimuli[count1:finalN,]
            CNOdiv$valueInhibitors = CNOdiv$valueInhibitors[count1:finalN,]
            for(b in 1:length(CNOdiv$valueSignals)) {
                CNOdiv$valueSignals[[b]] = CNOdiv$valueSignals[[b]][count1:finalN,]
            }
            for(d in 1:length(simDiv)) {
                simDiv[[d]] = simDiv[[d]][count1:finalN,]
            }
            count1 = count1 + div1
            CNOlistSet = c(CNOlistSet, list(CNOdiv))
            simResultsSet = c(simResultsSet, list(simDiv))
        }
    } else {

        CNOlistSet = list(CNOlist)
        simResultsSet = list(simResults)
    }

    for(f in 1:length(CNOlistSet)) {

            plotOptimResultsPan(
            simResults=simResultsSet[[f]],
            CNOlist=CNOlistSet[[f]],
            formalism="ss2",
            tPt=tPt,plotParams=plotParams
            )

        if(plotPDF == TRUE) {
            if(is.null(tag)) {
                filename <- paste("SimResultsT1T2", f, ".pdf", sep="")
            } else {
                filename <- paste(tag,"SimResultsT1T2", f, ".pdf", sep="_")
            }
            plotOptimResultsPan(
                simResults=simResultsSet[[f]],
                CNOlist=CNOlistSet[[f]],
                pdf=TRUE,
                formalism="ss2",
                pdfFileName=filename,
                tPt=tPt, plotParams=plotParams
            )
        }
    }
}


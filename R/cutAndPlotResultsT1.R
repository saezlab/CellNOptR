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
# $Id: cutAndPlotResultsT1.R 2256 2012-08-29 15:49:30Z cokelaer $

cutAndPlotResultsT1 <- function(model, bString, simList=NULL, CNOlist, indexList=NULL,
 plotPDF=FALSE, tag=NULL, tPt=CNOlist$timeSignals[2], plotParams=list(maxrow=10))
{

    if (is.null(simList)==TRUE){
        simList = prep4sim(model)
    }
    if (is.null(indexList)==TRUE){
        indexList = indexFinder(CNOlist, model)
    }

    if ("maxrow" %in% names(plotParams) == FALSE){
        plotParams$maxrow = 10
    }


    # keep simList and indxList for back compatibility ?
    modelCut <- cutModel(model, bString)
    simListCut <- cutSimList(simList, bString)

    # t0
    Sim0 <- simulatorT0(CNOlist=CNOlist, model=modelCut, simList=simListCut, indexList=indexList)
    simRes0 <- as.matrix(Sim0[,indexList$signals])

    # t1
    Sim <- simulatorT1(CNOlist=CNOlist, model=modelCut, simList=simListCut, indexList=indexList)
    simRes <- as.matrix(Sim[,indexList$signals])

    simResults <- list(t0=simRes0, t1=simRes)

    # if there is a lot of data, split up cnolist
    # make the max dimensions 10 x 10

    dim1 = dim(CNOlist$valueSignals[[1]])[1]
    dim2 = dim(CNOlist$valueSignals[[1]])[2]

    CNOlistSet = list()
    simResultsSet = list()

    if(dim1 > plotParams$maxrow) { #|| dim2 > 10) {

        par1 = ceiling(dim1/plotParams$maxrow)
        div1 = ceiling(dim1/par1)
        #par2 = ceiling(dim2/plotParams$maxrow)
        #div2 = ceiling(dim2/par2)

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

    outputFilenames = list()
    for(f in 1:length(CNOlistSet)) {

        plotOptimResultsPan(
            simResults=simResultsSet[[f]],
            CNOlist=CNOlistSet[[f]],
            formalism="ss1",
            tPt=tPt,
            plotParams=plotParams
            )

        if(plotPDF == TRUE) {
            if(is.null(tag)) {
                filename <- paste("SimResultsT1_", f, ".pdf", sep="")
            } else {
                filename <- paste(tag,"SimResultsT1",f,".pdf",sep="_")
            }
            plotOptimResultsPan(
                simResults=simResultsSet[[f]],
                CNOlist=CNOlistSet[[f]],
                pdf=TRUE,
                formalism="ss1",
                pdfFileName=filename,
                tPt=tPt,
                plotParams=plotParams
            )
            outputFilenames[[f]] = filename
        }
    }
    return(outputFilenames)
}


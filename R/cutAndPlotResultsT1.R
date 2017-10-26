#
#  This file is part of the CNO software
#
#  Copyright (c) 2011-2012 - EMBL - European Bioinformatics Institute
#
#  File author(s): CNO developers (cno-dev@ebi.ac.uk)
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  CNO website: http://www.cellnopt.org
#
##############################################################################
# $Id$

cutAndPlotResultsT1 <- function(model, bString, simList=NULL, CNOlist, indexList=NULL,
 plotPDF=FALSE, tag=NULL, tPt=CNOlist@timepoints[2], plotParams=list(maxrow=10))
{

    if ((class(CNOlist)=="CNOlist")==FALSE){
        CNOlist = CellNOptR::CNOlist(CNOlist)
    } 
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
    simRes0 <- as.matrix(Sim0[,indexList$signals,drop=F])
    #simRes0 = Sim0
    # t1
    Sim <- simulatorT1(CNOlist=CNOlist, model=modelCut, simList=simListCut, indexList=indexList)
    simRes <- as.matrix(Sim[,indexList$signals,drop=F])
    #simRes = Sim

    simResults <- list(t0=simRes0, t1=simRes)

    # if there is a lot of data, split up cnolist
    # make the max dimensions 10 x 10

    dim1 = dim(CNOlist@signals[[1]])[1]
    dim2 = dim(CNOlist@signals[[1]])[2]

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
            CNOdiv@cues = CNOdiv@cues[count1:finalN,,drop=F]
            CNOdiv@stimuli = CNOdiv@stimuli[count1:finalN,,drop=F]
            CNOdiv@inhibitors = CNOdiv@inhibitors[count1:finalN,,drop=F]
            for(b in 1:length(CNOdiv@signals)) {
                CNOdiv@signals[[b]] = CNOdiv@signals[[b]][count1:finalN,,drop=F]
            }
            for(d in 1:length(simDiv)) {
                simDiv[[d]] = simDiv[[d]][count1:finalN,,drop=F]
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

        mse = plotOptimResultsPan(
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
            mse = plotOptimResultsPan(
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
    return(list(filenames=outputFilenames, mse=mse, simResults=simResultsSet))
}


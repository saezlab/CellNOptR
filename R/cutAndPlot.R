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

cutAndPlot <- function(CNOlist, model, bStrings, plotPDF=FALSE, tag=NULL,
 plotParams=list(maxrow=10))
{

    # bitStrings must be a list of bitString (T1, T2, ...TN)
    tPt = length(bStrings)+1

    simList <- prep4sim(model)
    indexList <- indexFinder(CNOlist=CNOlist,model=model)

    

    # if tPt nothing to plot
    if (tPt==1){
        stop("noting to do with time data at time 0")
    }

    # if tPt=2 (default), call cutAndPlotResultsT1
    if (tPt == 2){
       #print("Entering cutAndPlotResultsT1")

       outputs = cutAndPlotResultsT1(model=model, bString=bStrings[[1]], simList=simList, 
            CNOlist=CNOlist, indexList=indexList, plotPDF=plotPDF, tag=tag,
            plotParams=plotParams)
    }

    if (tPt>=3){
       #print("Entering cutAndPlotResultsTN")
       outputs = cutAndPlotResultsTN(
         CNOlist=CNOlist,
         model=model,
         bStrings=bStrings,
         plotPDF=plotPDF,
         tag=tag) 
    }
    # if tPt=2 (default), call cutAndPlotResultsT2
    #return(outputs)
}

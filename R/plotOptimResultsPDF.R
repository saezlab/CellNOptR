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
# $Id: plotOptimResultsPDF.R 3155 2013-01-09 15:24:58Z cokelaer $
plotOptimResultsPDF<-function(
    simResults=simResults,
    expResults=expResults,
    times=times,
    namesCues=namesCues,
    namesSignals=namesSignals,
    valueCues=valueCues,
    filename, formalism="new"){

    if(sum(dim(simResults[[1]])) < 20){

        pdf(file=filename,width=14,height=7)

        }else{

            pdf(file=filename,width=21,height=10)

            }

    plotOptimResults(
        simResults=simResults,
        expResults=expResults,
        times=times,
        namesCues=namesCues,
        namesSignals=namesSignals,
        valueCues=valueCues, 
        formalism=formalism)

    dev.off()
    }


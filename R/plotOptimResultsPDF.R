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
# $Id: plotOptimResultsPDF.R 2243 2012-08-28 15:35:11Z cokelaer $
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


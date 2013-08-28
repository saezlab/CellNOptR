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
# $Id: plotFit.R 3522 2013-04-23 11:01:03Z cokelaer $
plotFit<-function(optRes, filename=NULL){

   oldPar = par(no.readonly = TRUE)
   on.exit(par(oldPar))

   if (is.null(filename)!=TRUE){
        pdf(filename)
    }

    par(mfrow=c(2,1),mar=c(0.5,4,4,0))
    plot(
        x=optRes$results[,"Generation"],
        y=optRes$results[,"Avg_Score_Gen"],
        xlab=NA,
        xaxt="n",
        ylab="Average score of generation",
        type="l")
    par(mar=c(4,4,0,0))
    plot(
        x=optRes$results[,"Generation"],
        y=optRes$results[,"Best_score_Gen"],
        xlab="Generations",
        ylab="Best Score",
        type="l")

    if (is.null(filename)!=TRUE){
        dev.off()
    }

   # par(oldPar)

} 


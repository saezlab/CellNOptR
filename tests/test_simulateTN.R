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
# This test runs T3 data set to compute the scores (must be zero), and run
# cutAndPlot plot on each time points.

# DO NOT MODIFY please (TC, June 2012)

library(CellNOptR)
pknmodel = readSIF(system.file("ToyModelT3/ToyModelT3.sif", package="CellNOptR"))
cnolist = CNOlist(system.file("ToyModelT3/ToyDataT3.csv",  package="CellNOptR"))

model = preprocessing(cnolist, pknmodel)

# computeScoreT1 with init string made of ones
verbose = FALSE
initBstring<-rep(1,length(model$reacID))
score = computeScoreT1(cnolist, model,bString=rep(1,length(model$reacID)))
bestBS = c(1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 0,0)
bestBS2 <- c(0,0,0,0,0,0,1) # given by running gaBinaryTN
bestBS3 <- c(0,0,1,0,0,0) # given by running gaBinaryTN


SimT1<-simulateTN(CNOlist=cnolist, model=model, bString=list(bestBS))
SimT2<-simulateTN(CNOlist=cnolist, model=model, bStrings=list(bestBS, bestBS2))
SimT3<-simulateTN(CNOlist=cnolist, model=model, bStrings=list(bestBS, bestBS2, bestBS3))

score1 = computeScoreT1(cnolist, model, bString=bestBS)
score2 = computeScoreTN(cnolist, model, bStrings=list(bestBS,bestBS2))
score3 = computeScoreTN(cnolist, model, bStrings=list(bestBS,bestBS2, bestBS3))

print(score1)
print(score2)
print(score3)
if (score1>0.01 || score2>0.024){
   stop("errore")
}
cutAndPlot(cnolist, model,bStrings=list(bestBS),plotPDF=TRUE, tag="test1")
cutAndPlot(cnolist, model,bStrings=list(bestBS,bestBS2),plotPDF=TRUE, tag="test2")
cutAndPlot(cnolist, model, bStrings=list(bestBS,bestBS2,bestBS3),plotPDF=TRUE, tag="test3")


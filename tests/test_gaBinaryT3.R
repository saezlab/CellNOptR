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
#  CNO website: http://www.ebi.ac.uk/saezrodriguez/cno
#
##############################################################################
# $Id$

# This is a test of the ToyModel and gaBinaryT1
# 

library(CellNOptR)

pknmodel = readSIF(system.file("ToyModelT3/ToyModelT3.sif", package="CellNOptR"))
cnolist = CNOlist(system.file("ToyModelT3/ToyDataT3.csv", package="CellNOptR"))
model = preprocessing(cnolist, pknmodel, verbose=FALSE)

# expected values
bestBS = c(1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 0,0)
bestBS2 <- c(0,0,0,0,0,0,1)
bestBS3 <- c(0,0,1,0,0,0)

# just to check that the simulateTN function works 
SimT1<-simulateTN(CNOlist=cnolist, model=model, bString=list(bestBS))
SimT2<-simulateTN(CNOlist=cnolist, model=model, bStrings=list(bestBS,bestBS2))
SimT3<-simulateTN(CNOlist=cnolist, model=model, bStrings=list(bestBS, bestBS2, bestBS3))


# again, just to check that gaBinary works
# run T1 first, 
T1opt<-gaBinaryT1(CNOlist=cnolist,model=model,verbose=FALSE)
# run T2
T2opt<-gaBinaryTN(CNOlist=cnolist,model=model,bStrings=list(T1opt$bString),verbose=FALSE)
# run T3
T3opt<-gaBinaryTN(CNOlist=cnolist,model=model,bStrings=list(bestBS, bestBS2),verbose=FALSE)

print( T1opt$bScore)
print( T2opt$bScore)
print( T3opt$bScore)
# no using the hardcoded parameters, we can check the output of te scores that
# must be tiny.
score1 = computeScoreT1(cnolist, model, bString=bestBS)
score2 = computeScoreTN(cnolist, model, bStrings=list(bestBS,bestBS2))
score3 = computeScoreTN(cnolist, model, bStrings=list(bestBS,bestBS2, bestBS3))


print(score1)
print(score2)
print(score3)
if (score1>0.01 || score2>0.24 || score3>0.1){
   # ideally, the score should all be close to 0. In practice, it's about 1e-5
   # However, in the  T3 case, once in while, the score is 0.0953 hence the
   # score3>0.1
   stop("errore")
}
cutAndPlot(cnolist, model,bStrings=list(bestBS),plotPDF=TRUE, tag="test1")
cutAndPlot(cnolist, model,bStrings=list(bestBS,bestBS2),plotPDF=TRUE, tag="test2")
cutAndPlot(cnolist, model, bStrings=list(bestBS,bestBS2,bestBS3),plotPDF=TRUE, tag="test3")



warnings()

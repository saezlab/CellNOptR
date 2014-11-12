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
# test if the priorBitString is used correctly in gaBinaryT1DO NOT MODIFY please (TC, June 2012)

library(CellNOptR)


# Load T1 and T2 data. Only T1 is used here below.
data(ToyModel, package="CellNOptR")
data(CNOlistToy, package="CellNOptR")
cnolist = CNOlistToy
pknmodel = ToyModel

# preprocessing
model = preprocessing(cnolist, pknmodel)

# computeScoreT1 with init string made of ones
initBstring<-rep(1,length(model$reacID))



priorBitString = rep(NA, length(model$reacID))
priorBitString[1] = 0
priorBitString[2] = 0
priorBitString[3] = 0
priorBitString[4] = 0

# Second, you call the gaBinaryT1 function by providing the priorBitString
# argument:
ToyT1opt<-gaBinaryT1(CNOlist=CNOlistToy, model=model,
    initBstring=initBstring, maxGens=10, popSize=5,
    verbose=FALSE, priorBitString=priorBitString)


for (x in ToyT1opt$results[,7]){
     x = strsplit(x ,",")[[1]]
     if (as.numeric(x[[1]]) != 0 ){ stop("first element must be 0")}
     if (as.numeric(x[[2]]) != 0 ){ stop("second element must be 0")}
     if (as.numeric(x[[3]]) != 0 ){ stop("tirdh element must be 0")}
     if (as.numeric(x[[4]]) != 0 ){ stop("fourth element must be 0")}
}



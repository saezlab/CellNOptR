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
library(CellNOptR)

# read data
pknmodel = readSIF(system.file("ToyModel/ToyPKNMMB.sif", package="CellNOptR"))
data = readMIDAS(system.file("ToyModel/ToyDataMMB.csv", package="CellNOptR"))
cnolist = makeCNOlist(data, subfield=FALSE)

# preprocessing
model = preprocessing(cnolist, pknmodel, verbose=FALSE)

# optimisation
T1opt<-gaBinaryT1(CNOlist=cnolist,model=model,verbose=FALSE)

truebs =c(1, 1, 0, 1, 1, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0)
print(truebs)
print(T1opt$bString)

# testing valid output
if (dist(rbind(T1opt$bString, truebs))>2){
stop("something wrong going on")
}

# extra call to simulateTN
SimT1<-simulateTN(CNOlist=cnolist,model=model, bStrings=list(truebs))


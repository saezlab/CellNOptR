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
data = readMIDAS(system.file("ToyModelT3/ToyDataT3.csv", package="CellNOptR"))
cnolist = makeCNOlist(data, subfield=FALSE)
model = preprocessing(cnolist, pknmodel, verbose=FALSE)

# expected values
truebs = c(1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 0,0)
truebs2 <- c(0,0,0,0,0,0,1)
truebs3 <- c(0,0,1,0,0,0)
truebs3_bis <- c(0,0,0,1,0,0)
truebs3_ter <- c(0,0,0,0,1,0)


# run T1 first, 
T1opt<-gaBinaryT1(CNOlist=cnolist,model=model,verbose=FALSE)
print(T1opt$bString)


if (all(T1opt$bString == truebs)==FALSE){
stop("something wrong going on")
}


# run T2
T2opt<-gaBinaryTN(CNOlist=cnolist,model=model,bStrings=list(T1opt$bString),verbose=FALSE)
print(T2opt$bString)
if (all(T2opt$bString == truebs2)==FALSE){
    stop("something wrong going on T2")
}

# run T3
T3opt<-gaBinaryTN(CNOlist=cnolist,model=model,bStrings=list(truebs, truebs2),verbose=FALSE)
print(T3opt$bString)


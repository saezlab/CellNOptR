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
data(CNOlistToy2, package="CellNOptR")
data(ToyModel2, package="CellNOptR")
cnolist = CNOlistToy2
pknmodel = ToyModel2

model = preprocessing(cnolist, pknmodel, verbose=FALSE)
T1opt<-gaBinaryT1(CNOlist=cnolist,model=model,verbose=FALSE)
print(T1opt$bString)

truebs = c(1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0)
if (all(T1opt$bString == truebs)==FALSE){
    stop("something wrong going on")
}


T2opt<-gaBinaryTN(CNOlist=cnolist,model=model,bStrings=list(T1opt$bString),verbose=FALSE)
print(T2opt$bString)
truebs2 = c(0, 0, 0, 1, 0, 0, 0)
if (all(T2opt$bString == truebs2)==FALSE){
    stop("something wrong going on")
}


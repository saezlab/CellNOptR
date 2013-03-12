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

# This test runs optimisation using Boolean logic at different time points
# Author: S.Schrier, T.Cokelaer



manySteadyStates <-function(
  CNOlist,
  model,
  sizeFac=0.0001,
  NAFac=1,
  popSize=50,
  pMutation=0.5,
  maxTime=60,
  maxGens=500,
  stallGenMax=100,
  selPress=1.2,
  elitism=5,
  relTol=0.1,
  verbose=FALSE,
  priorBitString=NULL){


    #initialize Opt array
    Opt<-list()    
    #initialze a bString array  
    bStrings<-list()
    simRes<-list()


    T1opt<-gaBinaryT1(CNOlist=CNOlist,
                  model=model,
                  stallGenMax=stallGenMax,
                  sizeFac=sizeFac,
                  NAFac=NAFac,
                  popSize=popSize,
                  pMutation=pMutation,
                  maxTime=maxTime,
                  maxGens=maxGens,
                  selPress=selPress,
                  elitism=elitism,
                  relTol=relTol,
                  verbose=verbose,
                  priorBitString=priorBitString)

    Opt[[1]]<-T1opt
    simT1<-simulateTN(CNOlist=CNOlist, model=model, bStrings=list(T1opt$bString))
    simRes[[1]]<-simT1
    bStrings[[1]] = T1opt$bString

    if (length(CNOlist@signals)>2){
        for(i in 3:length(CNOlist@signals)){
            Opt[[i-1]]<-gaBinaryTN(CNOlist=CNOlist,
                        model=model,
                        bStrings=bStrings,
                        stallGenMax=stallGenMax,
                        maxTime=maxTime,
                        sizeFac=sizeFac,
                        NAFac=NAFac,
                        popSize=popSize,
                        pMutation=pMutation,
                        maxGens=maxGens,
                        selPress=selPress,
                        elitism=elitism,
                        relTol=relTol,
                        verbose=verbose,
                        priorBitString=priorBitString)

            bStrings[[i-1]] = Opt[[i-1]]$bString

            simRes[[i]]<-simulateTN(CNOlist,model,bStrings)
        }
    }   
    return(list(bStrings=bStrings, Opt=Opt, simRes=simRes))
}



library(CellNOptR)

# one steady state
data(CNOlistToy, package="CellNOptR")
cnolist = CNOlist(CNOlistToy)
data(ToyModel, package="CellNOptR")
results = manySteadyStates(cnolist, ToyModel)
print(results$bStrings)

# two steady state
data(CNOlistToy2, package="CellNOptR")
cnolist = CNOlist(CNOlistToy2)
data(ToyModel2, package="CellNOptR")
results = manySteadyStates(cnolist, ToyModel2)
print(results$bStrings)

# 3 steady state
ToyModel3 = readSIF(system.file("ToyModelT3/ToyModelT3.sif", package="CellNOptR"))
cnolist = CNOlist(system.file("ToyModelT3/ToyDataT3.csv", package="CellNOptR"))
results = manySteadyStates(cnolist, ToyModel3)
print(results$bStrings)



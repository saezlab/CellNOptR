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

CNORbool<-function(CNOlist, model, paramsList=defaultParameters(),
    compression=TRUE, expansion=TRUE, cutNONC=TRUE, verbose=FALSE)
{
    # just an alis
    cnolist = CNOlist

    # if filenames provided, need to call the proper functions
    if (is.character(cnolist)==TRUE){
        cnolist = makeCNOlist(readMIDAS(cnolist), subfield=FALSE)
    }
    if (is.character(model)==TRUE){
        model = readSIF(model)
    }

    # overwrite verbosity
    paramsList$verbose = verbose

    # 1. Checks data to model compatibility
    checkSignals(cnolist, model)
    # 2.preprocessing
    model = preprocessing(cnolist, model, compression=compression,
        expansion=expansion, cutNONC=cutNONC, verbose=verbose)

    #8.Optimisation t1
    bStrings = list()
    initBstring<-rep(1,length(model$reacID))
    print("Entering gaBinaryT1")
    T1opt<-gaBinaryT1(CNOlist=cnolist,
        model=model,
        initBstring=initBstring,
        sizeFac=paramsList$sizeFac,
        NAFac=paramsList$NAFac,
        popSize=paramsList$popSize,
        pMutation=paramsList$pMutation,
        maxTime=paramsList$maxTime,
        maxGens=paramsList$maxGens,
        stallGenMax=paramsList$stallGenMax,
        selPress=paramsList$selPress,
        elitism=paramsList$elitism,
        relTol=paramsList$relTol,
        verbose=paramsList$verbose)
    bStrings[[1]] = T1opt$bString


    if (length(cnolist$valueSignals)==2){
        return(list(model=model, bStrings=bStrings))
    }

    #.Optimise tN where N>1
    Times = 1
    T2opt<-NA # default value

    for (i in 3:length(cnolist$valueSignals)){
        Times = Times + 1
        print(paste("Entering gaBinaryTN (time=", Times, ")", sep=" "))

        TNopt<-gaBinaryTN(
            CNOlist=cnolist,
            model=model,
            bStrings=bStrings,
            sizeFac=paramsList$sizeFac,
            NAFac=paramsList$NAFac,
            popSize=paramsList$popSize,
            pMutation=paramsList$pMutation,
            maxTime=paramsList$maxTime,
            maxGens=paramsList$maxGens,
            stallGenMax=paramsList$stallGenMax,
            selPress=paramsList$selPress,
            elitism=paramsList$elitism,
            relTol=paramsList$relTol,
            verbose=paramsList$verbose)

        bStrings[[Times]] = TNopt$bString
        print(TNopt$bString)

        if (Times==2){
            T2opt = TNopt
        }
    }

    return(list(model=model, bStrings=bStrings))
}

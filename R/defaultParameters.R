#
#  This file is part of the CNO software
#
#  Copyright (c) 2011-2012 - EMBL - European Bioinformatics Institute - Massachusetts Institute of Technology
#
#  File author(s): CNO developers (cno-dev@ebi.ac.uk)
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  CNO website: http://www.ebi.ac.uk/saezrodriguez/software/cno
#
##############################################################################
# $Id$

defaultParameters<-function(data=NA, model=NA){

    paramsList<-list()
    paramsList$data<-data
    paramsList$model<-model

    # GA parameters
    paramsList$sizeFac<-1e-04
    paramsList$NAFac<-1
    paramsList$popSize<-50
    paramsList$pMutation<-0.5
    paramsList$maxTime<-3*60
    paramsList$maxGens<-500
    paramsList$stallGenMax<-100
    paramsList$selPress<-1.2
    paramsList$elitism<-5
    paramsList$relTol<-0.1
    paramsList$verbose<-FALSE

    return(paramsList)
}


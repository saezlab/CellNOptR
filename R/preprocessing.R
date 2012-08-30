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

preprocessing<-function(data=NULL, model, cutNONC=TRUE, compression=TRUE,
    expansion=TRUE, ignoreList=NA, maxInputsPerGate=2,verbose=TRUE){

    # why not doing this check here ? Does not cost too much
    if (is.null(data)!=TRUE){
        checkSignals(CNOlist=data,model=model)
    }

    # a copy of the model
    cutModel <- model

    if (cutNONC==TRUE && is.null(data)!=TRUE){
        # Find the indices, in the model, of the species that are inh/stim/sign
        indices<-indexFinder(CNOlist=data, model=model,    verbose=verbose)

        # Find the indices of the non-osb/non-contr
        temp_indices <- findNONC(model=model, indexes=indices,verbose=verbose)
        # Cut the nonc off the model
        cutModel <-cutNONC(model=model, NONCindexes=temp_indices)
    }

    if (compression == TRUE && is.null(data)!=TRUE){
        # Recompute the indices
        temp_indices<-indexFinder(CNOlist=data, model=cutModel)

        # Compress the model
        cutModel<-compressModel(model=cutModel,indexes=temp_indices)
    }

    # Recompute the indices. We can do it now because the expanson gate does not
    # remove species but only add and gates.
    #if (is.null(data)!=TRUE){
    #    indices<-indexFinder(CNOlist=data,model=cutModel)
    #}
    #else{
    #    indices <- NULL
    #}

    # Expand the gates
    if (expansion == TRUE){
        cutModel <- expandGates(model=cutModel, ignoreList=ignoreList,maxInputsPerGate=maxInputsPerGate)
    }

    # since version 1.3.28 return only model, indices are recomputed in other
    # functions
    return(cutModel)
}

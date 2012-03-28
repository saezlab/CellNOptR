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

preprocessing<-function(Data, Model, cutnonc=TRUE, compression=TRUE,
expansion=TRUE, ignoreList=NA, maxInputsPerGate=2,verbose=TRUE){

    # why not doing this check here ? Does not cost too much
	checkSignals(CNOlist=Data,Model=Model)

    # a copy of the model 
    cutmodel <- Model

    if (cutnonc==TRUE){
        # Find the indices, in the model, of the species that are inh/stim/sign
	    indices<-indexFinder(CNOlist=Data, Model=Model,	verbose=verbose)

        # Find the indices of the non-osb/non-contr	
	    temp_indices <- findNONC(Model=Model, indexes=indices,verbose=verbose)
        # Cut the nonc off the model
        cutmodel <-cutNONC(Model=Model, NONCindexes=temp_indices)
    }

    if (compression == TRUE){
        # Recompute the indices
	    temp_indices<-indexFinder(CNOlist=Data, Model=cutmodel)

        # Compress the model
    	cutmodel<-compressModel(Model=cutmodel,indexes=temp_indices)
    }

    # Recompute the indices. We can do it now because the expanson gate does not
    # remove species but only add and gates.
	indices<-indexFinder(CNOlist=Data,Model=cutmodel)

    # Expand the gates	
    if (expansion == TRUE){
        cutmodel <- expandGates(Model=cutmodel, ignoreList=ignoreList,maxInputsPerGate=maxInputsPerGate)
    }

    return( list(model=cutmodel, indices=indices))
}

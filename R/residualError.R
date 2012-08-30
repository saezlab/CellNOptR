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
# $Id: residualError.R 2267 2012-08-30 15:31:54Z cokelaer $
residualError<-function(CNOlist){


    if ((class(CNOlist)=="CNOlist")==FALSE){
         CNOlist = CellNOptR::CNOlist(CNOlist)
     }

    nTimes = length(CNOlist@timepoints) - 1  # we do not take into account t0


	resErr<-rep(NA, nTimes + 2) # for back compatibility, we will store t1andt2
                                # (+1) but we also store the total (+1)

    # build up the name of the columns: t1, t2, ...tN, t1andt2, total
    namesresErr = NULL
    for (i in 1:nTimes){
        namesresErr = cbind(namesresErr, paste("t", i, sep=""))
    }
    namesresErr = cbind(namesresErr, "t1andt2")
    namesresErr = cbind(namesresErr, "total")
    names(resErr) <- namesresErr

    # compute errors at each time 
    for (i in 1:nTimes){
        # we skip time 0. t1 starts at index 2 hence the i+1
        Diff <- round(CNOlist@signals[[i+1]])-CNOlist@signals[[i+1]]
        resErr[i]<-sum(Diff^2, na.rm=TRUE)
    }

    # compute the total and t1andt2 columns
    resErr["total"]<-sum(resErr, na.rm=T) # should be before t1ant2 being filled
	if(length(CNOlist@signals) >= 3){
        resErr[nTimes+1]<-resErr["t1"]+resErr["t2"]
    }

	return(resErr)	
	}


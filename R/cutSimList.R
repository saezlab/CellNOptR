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
cutSimList <- function(SimList, bitString)
{
    bitString = as.logical(bitString)
    SimListCut <- SimList
    finalCube <- SimListCut$finalCube[bitString,]
    ixNeg <-SimListCut$ixNeg[bitString,]
    ignoreCube <- SimListCut$ignoreCube[bitString,]
    maxIx <- SimListCut$maxIx[bitString]
    # in some cases the finalcube is a matrix but list of integer, so we
    # need to convert back to a matrix. Happens for simple models only hence
    # the warning.
    if (is.matrix(finalCube) == FALSE){
        #warning("converting back to matrix in prep4sim")
        SimListCut$finalCube<-matrix(finalCube,
            dimnames=list(names(finalCube), 1))
        SimListCut$ixNeg<-matrix(ixNeg, dimnames=list(names(ixNeg), 1))
        SimListCut$ignoreCub<-matrix(ignoreCube,dimnames=list(names(ignoreCube), 1))
        SimListCut$maxIx<-matrix(maxIx,dimnames=list(names(maxIx), 1))
    }
    else{
        SimListCut$finalCube = finalCube
        SimListCut$ixNeg<-ixNeg
        SimListCut$ignoreCube<-ignoreCube
        SimListCut$maxIx<-maxIx
    }
    return(SimListCut)
}

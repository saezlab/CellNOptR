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


writeSIF <- function(model, filename){

    # convert internal model structure to a SIF matrix
    sif = model2sif(model)

    # and save it into a file
    if (file.exists(filename)==FALSE){
        write.table(sif[,1:3],file=filename,
            row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
    }
    else{
       stop(paste("File ", filename, "already exists.",  sep=""))
    }


}


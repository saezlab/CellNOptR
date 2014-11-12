#
#  This file is part of the CNO software
#
#  Copyright (c) 2011-2012 - EMBL - European Bioinformatics Institute
#
#  File author(s): CNO developers (cno-dev@ebi.ac.uk)
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  CNO website: http://www.cellnopt.org
#
##############################################################################
# $Id$


writeSIF <- function(model, filename, overwrite=FALSE){

    # convert internal model structure to a SIF matrix
    sif = model2sif(model)

    # and save it into a file
    if (file.exists(filename)==FALSE){
        write.table(sif[,1:3],file=filename,
            row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
    }
    else{
       if (overwrite==FALSE){
            write.table(sif[,1:3],file=filename,
                row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
        }
        else{
           stop(paste("File ", filename, "already exists.",  sep=""))
        }
    }


}


#
#  This file is part of the CNO software
#
#  Copyright (c) 2011-2012 - EBI
#
#  File author(s): CNO developers (cno-dev@ebi.ac.uk)
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  CNO website: http://www.ebi.ac.uk/saezrodriguez/cno
#
##############################################################################
# $Id$
CNOdata <- function(filename, verbose=FALSE, url=NULL){


    valid_filenames = c(
        "PKN-ToyMMB.sif", "PKN-ToyPB.sif", "PKN-ToyMSB2009.sif",
        "PKN-ToyMMB_T2.sif", "MD-ToyMMB_T2.csv",
        "PKN-ExtLiverPCB.sif", "MD-ToyMMB.csv", "MD-ToyPB.csv", 
        "MD-ToyMSB2009.csv", "MD-ExtLiverPCB.csv", "PKN-ToyPCB.sif",
         "MD-ToyPCB.csv")


    if ((filename %in% valid_filenames)==FALSE){
        print("Provided filename not registered. Please use one of ")
        print(valid_filenames)
        stop()
    }

    #directory = strsplit(filename, split="-", fixed=TRUE)[[1]][2]
    #print(directory)
    #directory = strsplit(directory, split=".", fixed=TRUE)[[1]][1]
    #print(directory)

    if (is.null(url)==TRUE){
        url = "http://www.ebi.ac.uk/~cokelaer/cellnopt/data/_downloads"
    }

    #filename = paste(URL, directory, filename, sep="/")
    filename = paste(url,  filename, sep="/")
    if (verbose==TRUE){
        print(filename)
    }
    library(RCurl)
    data = getURL(filename)
    fh = tempfile("cellnopt_", fileext=".dat")
    if (verbose==TRUE){
        print(paste("data downloaded and copied into ", fh, sep=" "))
    }
    write(data, fh)

    return (fh)

}

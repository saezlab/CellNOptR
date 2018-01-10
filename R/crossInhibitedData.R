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
#  CNO website: http://www.ebi.ac.uk/saezrodriguez/cno
#
##############################################################################
# $Id$

crossInhibitedData <- function(object){
          times = object@timepoints 

          # identify names found in inhibitors and signals list
          inhibitors = colnames(cnolist@inhibitors)[colnames(cnolist@inhibitors) %in% colnames(cnolist@signals[[1]])]

          # only those ones must be crossed
          for (inhibitor in inhibitors){
              for (time in seq_along(times)){
                  mask = cnolist@inhibitors[,inhibitor] == 1
                  object@signals[[time]][mask, inhibitor] = NA
              }
          }
    return(object)
}



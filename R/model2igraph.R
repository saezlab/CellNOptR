#
#  This file is part of the CNO software
#
#  Copyright (c) 2011-2013 - EBI
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

model2igraph <-function(model){

    # The input is a model returned by readSIF()
    # need to convert it back to a sif and then a graphNEL object
    g = sif2graph(model2sif(model))

    # now we can get a igraph object
    library(igraph)
    g2 = igraph.from.graphNEL(g)

    # so it can be saved in many different format.
    return(g2)
}

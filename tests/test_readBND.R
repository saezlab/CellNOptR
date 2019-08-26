#
#  This file is part of the CNO software
#
#  Copyright (c) 2011-2019 - EBI
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
# testing BND imports

library(CellNOptR)

bnd_network = readBND(system.file("test/example.bnd", package="CellNOptR"))

if(!all(c("reacID", "namesSpecies", "interMat", "notMat" ) %in% names(bnd_network))){
	stop("missing field")
}

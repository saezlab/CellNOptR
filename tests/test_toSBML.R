
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
# test SBML export


library(CellNOptR)

pknmodel = readSIF(system.file("ToyModelT3/ToyModelT3.sif", package="CellNOptR"))

toSBML(pknmodel,tempfile(pattern = "test_SBMLqual.xml"))


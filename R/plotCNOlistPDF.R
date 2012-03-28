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
# $Id: plotCNOlistPDF.R 491 2012-02-02 17:59:17Z cokelaer $
plotCNOlistPDF <-
function(CNOlist,filename){
	pdf(file=filename,width=14,height=7)
	plotCNOlist(CNOlist)
	dev.off()
	}


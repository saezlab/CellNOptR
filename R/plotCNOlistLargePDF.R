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
# $Id: plotCNOlistLargePDF.R 3155 2013-01-09 15:24:58Z cokelaer $
plotCNOlistLargePDF <-
function(CNOlist,filename,nsplit, width=14, height=7){
    pdf(file=filename,width=width,height=height)
    plotCNOlistLarge(CNOlist,nsplit=nsplit,newDevice=FALSE)
    dev.off()
    }


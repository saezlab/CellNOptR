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
plotCNOlist2<-function(cnolist, simulated_cnolist=NULL, ymin=0, ymax=1){

    #check that CNOlist is a CNOlist


    if (is.null(simulated_cnolist)==FALSE){
        d1 = .cnolist2ggdata(cnolist)
        d2 = .cnolist2ggdata(simulated_cnolist, label="simulated")
        d = rbind(d1, d2)
    } else{
        d = .cnolist2ggdata(cnolist)
    }

    times=NULL # just to prevent warnings
    type=NULL # just to prevent warnings
    # in ggplot, you can remove colour=type to get only black lines
    # remove shape=type to get black-circle symbol
    sp <- ggplot(d, aes(x=times, y=values, group=type, shape=type, colour=type)) +
        geom_point(size=3, fill="white") +
        geom_line(size=1.) +
        xlab("") + 
        ylab("Conditions")
        #scale_shape_manual(name="data type", values=c(22,21))
        #scale_color_manual(name="data type", values=c("black","red"))

    sp + facet_grid(conditions~species) + theme_bw()+ 
        theme(
            strip.text.x=element_text(size=8,angle=0),
            strip.text.y = element_text(size=12, face="bold"),
            strip.background = element_rect(colour="red", fill="#CCCCFF")) +
         ylim(ymin, ymax)
    #+ggtitle("fff")

    # does not seem to work well if the following code is uncommented, ggplot
    # does nto show the plot even though the code does ot enter in the if
    # statement.
    # save in PDF or whateve accepted format by ggsave.
    #if (is.null(filename)==FALSE){
        #ggsave(filename, useDingbats=FALSE)
    #}
    #cannot return anything otherwise ggplot does not show the plot...
    #return(d)
}

#convert a cnolist signal data to a data structure suitable for ggplot
.cnolist2ggdata <- function(cnolist, label="measured"){
    allvalues = c()
    alltimes  = c()
    allspecies = c()
    allconditions = c()
    allcues = c()
    times = cnolist@timepoints
    nCues = length(cnolist@cues[1,])

    for(r in 1:dim(cnolist@signals[[1]])[1]){
        for(c in 1:dim(cnolist@signals[[1]])[2]){
            yVal = c()
            for (time in 1:length(times)){
                yVal<- c(yVal, as.numeric(cnolist@signals[[time]][r,c]))
            }
            specy = colnames(cnolist@signals[[1]])[c]
            species = rep(specy, length(yVal))
            cues = cnolist@cues[r,]
            conditions=rep(r, length(yVal))

            allcues = c(allcues, cues)
            allvalues = c(allvalues, yVal)
            allspecies = c(allspecies, species)
            alltimes = c(alltimes, times)
            allconditions = c(allconditions, conditions)
        }
    }
    n1 = length(allvalues)
    d1 = data.frame(values=allvalues, species=allspecies, times=alltimes,
        conditions=allconditions, type=rep(label, n1))

    ma = matrix(allcues, ncol=nCues)
    colnames(ma) = names(cnolist@cues[1,])
    d1 = cbind(d1, ma)

    return(d1)
}


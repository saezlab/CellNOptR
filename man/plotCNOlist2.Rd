\name{plotCNOlist2}
\alias{plotCNOlist2}
\title{
Another version of plotCNOlist that allows to plot 2 cnolist in the same layout
to compare them. This function uses ggplot2 library. It is recommended for small
data sets (about 15 species). 
}
\description{
This function plots the data in a CNOlist as a matrix of plots with a  row for
each condition and a column for each signal, Cues are simply represented by a
number.
}
\usage{
plotCNOlist2(cnolist, simulated_cnolist=NULL, ymin=0,ymax=1)
}

\arguments{
  \item{cnolist}{
a CNOlist
}
  \item{simulated_cnolist}{
another cnolist
}
\item{ymin}{
    Change the lower y-limit (default is 0)
}
\item{ymax}{
    Change the lower y-limit (default is 1)
}


}
\details{
This function can plot either a single CNOlist, or 2 on top of each other.
}
\author{
T. Cokelaer
}
\seealso{
plotCNOlistPDF, plotCNOlistLarge, plotCNOlistLargePDF, plotCNOlist
}
\examples{
# this data set is not an object so we need to convert it
data(CNOlistToy,package="CellNOptR")
cnolist = CNOlist(CNOlistToy)
plotCNOlist2(cnolist)
}

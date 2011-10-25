\name{plotCNOlist}
\alias{plotCNOlist}
\title{
Plot the data in a  CNOlist
}
\description{
This function plots the data in a CNOlist as a matrix of plots with a  row for each condition and a column for each signal, and an extra plot for each row that specifies which cues are present.
}
\usage{
plotCNOlist(CNOlist)
}
\arguments{
  \item{CNOlist}{
a CNOlist
}
}
\details{
This function can plot the normalised values or the un-normalised ones, it just needs a CNOlist.
}
\value{
This function just produces a plot on your graphics window
}
\author{
C. Terfve
}
\seealso{
plotCNOlistPDF, plotCNOlistLarge, plotCNOlistLargePDF
}
\examples{
data(CNOlistToy,package="CellNOptR")
plotCNOlist(CNOlistToy)
}
\name{plotCNOlistLarge}
\alias{plotCNOlistLarge}
\title{
Plot the data in a  CNOlist, for lists with many conditions.
}
\description{
This function plots the data in a CNOlist as a matrix of plots with a  row for each condition and a column for each signal, and an extra plot for each row that specifies which cues are present.
}
\usage{
plotCNOlistLarge(CNOlist,nsplit=4)
}
\arguments{
  \item{CNOlist}{
a CNOlist
}
  \item{nsplit}{
the number of splits in the condition dimension (one new plot window will be produced for each split, i.e. if you have 80 conditions and specify 4 splits you will get 4 plots with 20 conditions each).
}
}
\details{
This function can plot normalised values or the un-normalised ones, it just needs a CNOlist. This function makes plots of CNOlists that are more readable when many conditions are present in the data.  In addition to plotting the conditions divided into multiple plots, this function also plots the cues divided in two columns, one for inhibitors and one for stimuli.
}
\value{
This function just produces plots on your graphics window.
}
\author{
C. Terfve
}
\seealso{
plotCNOlist, plotCNOlistPDF, plotCNOlistLargePDF
}
\examples{
data(CNOlistDREAM,package="CellNOptR")
plotCNOlistLarge(CNOlistDREAM, nsplit=2)
}

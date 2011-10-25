\name{plotCNOlistLargePDF}
\alias{plotCNOlistLargePDF}
\title{
Plots a CNOlist into a pdf file, for lists with many conditions.
}
\description{
This function is a wrapper for plotCNOlistLarge, that plots the output directly in a pdf file.
}
\usage{
plotCNOlistLargePDF(CNOlist, fileName, nsplit)
}
\arguments{
  \item{CNOlist}{
a CNOlist
}
  \item{fileName}{
a name for your pdf file, eg. "plot.pdf"
}
  \item{nsplit}{
the number os splits along the condition dimension (see plotCNOlistLarge)
}
}
\details{
This function makes plots of CNOlists that are more readable when many conditions are present in the data.  In addition to plotting the conditions divided into multiple plots, this function also plots the cues divided in two columns, one for inhibitors and one for stimuli.
}
\value{
This function doesn't return anything, it just produces a pdf file with your plots, in your current working directory.
}
\author{
C. Terfve
}
\seealso{
plotCNOlistLarge, plotCNOlist, plotCNOlistPDF
}
\examples{
tmpdir<-tempdir()
setwd(tmpdir)
data(CNOlistDREAM,package="CellNOptR")
plotCNOlistLargePDF(CNOlistDREAM, fileName="dreamData.pdf",nsplit=2)
}

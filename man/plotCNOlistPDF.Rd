\name{plotCNOlistPDF}
\alias{plotCNOlistPDF}
\title{
Plots a CNOlist into a pdf file.
}
\description{
This function is a wrapper for plotCNOlist, that plots the output directly in a pdf file.
}
\usage{
plotCNOlistPDF(CNOlist, fileName)
}
\arguments{
  \item{CNOlist}{
a CNOlist
}
  \item{fileName}{
a name for your pdf file, eg. "plot.pdf"
}
}
\value{
This function doesn't return anything, it just produces a pdf file containing your plot, in your working directory.
}
\author{
C. Terfve
}
\seealso{
plotCNolist, plotCNOlistLarge, plotCNOlistLargePDF
}
\examples{
tmpdir<-tempdir()
setwd(tmpdir)
data(CNOlistToy,package="CellNOptR")
plotCNOlistPDF(CNOlist=CNOlistToy,fileName="ToyModelGraph.pdf")
}


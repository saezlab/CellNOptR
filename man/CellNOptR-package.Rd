\name{CellNOptR-package}
\alias{CellNOptR-package}
\alias{CellNOptR}
\docType{package}
\title{
R version of CellNOptR, boolean features
}
\description{
This package does optimisation of boolean logic networks of signalling pathways based on a previous knowledge network and a set of data collected upon perturbation of some of the nodes in the network.
}
\details{
\tabular{ll}{
Package: \tab CellNOptR\cr
Type: \tab Package\cr
Version: \tab 1.10.0\cr
Date: \tab 2014-03-13\cr
License: \tab GPLv3\cr
LazyLoad: \tab yes\cr
}
}
\author{
T. Cokelaer, A. MacNamara, F. Eduati, S. Schrier, C. Terfve

Maintainer: T. Cokelaer <cokelaer@ebi.ac.uk>
}
\references{
J. Saez-Rodriguez, L. G. Alexopoulos, J. Epperlein, R. Samaga, D. A. Lauffenburger, S. Klamt and P. K. Sorger. Discrete logic modeling as a means to link protein signaling networks with functional analysis of
mammalian signal transduction, Molecular Systems Biology, 5:331, 2009.
}

\keyword{ package }
\examples{
library(CellNOptR)

# quick 1 time point optimisation of a Prior Knowledge Network to MIDAS data. 
data(CNOlistToy,package="CellNOptR")
data(ToyModel,package="CellNOptR")

pknmodel = ToyModel
cnolist = CNOlist(CNOlistToy)
model = preprocessing(cnolist, pknmodel)
results = gaBinaryT1(cnolist, model, verbose=FALSE)
plotFit(results)
cutAndPlot(cnolist, model, list(results$bString))

# Same as above and HTML report
CNORwrap(name="Toy", 
    namesData=list(CNOlist="ToyData",model="ToyModel"),
    data=CNOlistToy, model=pknmodel)
}


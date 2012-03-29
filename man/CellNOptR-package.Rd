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
Version: \tab 1.2.0\cr
Date: \tab 2011-03-29\cr
License: \tab GPLv2\cr
LazyLoad: \tab yes\cr
}
}
\author{
C.Terfve, T.Cokelaer


Maintainer: C.Terfve <terfve@ebi.ac.uk>
}
\references{
J. Saez-Rodriguez, L. G. Alexopoulos, J. Epperlein, R. Samaga, D. A. Lauffenburger, S. Klamt and P. K. Sorger. Discrete logic modeling as a means to link protein signaling networks with functional analysis of
mammalian signal transduction, Molecular Systems Biology, 5:331, 2009.
}

\keyword{ package }
\examples{
library(CellNOptR)

tmpdir<-tempdir()
setwd(tmpdir)

data(CNOlistToy,package="CellNOptR")
data(ToyModel,package="CellNOptR")

#From here there are 2 versions: 1. If you want to set the parameters yourself
pList<-list(
	Data=CNOlistToy,
	Model=ToyModel,
	sizeFac = 1e-04, 
	NAFac = 1, 
	PopSize = 10, 
	Pmutation = 0.5, 
	MaxTime = 60, 
	maxGens = 5, 
	StallGenMax = 5, 
	SelPress = 1.2, 
	elitism = 5, 
	RelTol = 0.1,
	verbose=TRUE)

CNORwrap(
	paramsList=pList,
	Name="Toy",
	NamesData=list(CNOlist="ToyData",Model="ToyModel"),
	Data=NA,
	Model=NA)

#2. If you want to keep the default parameters
\dontrun{
CNORwrap(
	paramsList=NA,
	Name="Toy",
	NamesData=list(CNOlist="ToyData",Model="ToyModel"),
	Data=CNOlistToy,
	Model=ToyModel)
}
}

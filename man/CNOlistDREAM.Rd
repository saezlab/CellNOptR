\name{CNOlistDREAM}
\docType{data}
\alias{CNOlistDREAM}
\title{Data used for the DREAM3 challenge}
\description{This data object contains the DREAM data used in the package vignette, already loaded and formatted as a CNOlist object. This is to be used with the model "DreamModel". This is a data collected on HepG2 cells cultivated with or without stimulation of tgfa, ilk, mek12, pi3k and p38, in combination with inhibition of igf1 and/or il1a.  Seven phosphoproteins are measured using Luminex xMAP assays: akt, erk12, ikb, jnk12, p38, hsp27 and mek12.}
\usage{CNOlistDREAM}
\format{CNOlistDREAM is a list with the fields "namesCues" (character vector), "namesStimuli" (character vector), "namesInhibitors" (character vector), "namesSignals" (character vector), "timeSignals" (numerical vector), "valueCues" (numerical matrix), "valueInhibitors" (numerical matrix), "valueStimuli"(numerical matrix), "valueSignals"(numerical matrix).}
\source{This data and model is extracted from the Matlab version of CellNOpt1.0 (\url{http://www.ebi.ac.uk/saezrodriguez/software.html#CellNetOptimizer}).}
\references{
\enumerate{
\item{J. Saez-Rodriguez, L. G. Alexopoulos, J. Epperlein, R. Samaga, D. A.
Lauffenburger, S. Klamt and P. K. Sorger. Discrete logic modeling as a means to
link protein signaling networks with functional analysis of mammalian signal
transduction, Molecular Systems Biology, 5:331, 2009.}
\item{ Prill RJ, Marbach D, Saez-Rodriguez J, Sorger PK, Alexopoulos LG, Xue X,
Clarke ND, Altan-Bonnet G, and Stolovitzky G. Towards a rigorous assessment of
systems biology models: the DREAM3 challenges. PLoS One, 5(2):e9202, 2010. }
}
}
\keyword{datasets}

\name{plotModel}
\alias{plotModel}
\title{
Plot a model
}
\description{
This function can be used to plot a prior model network before any
pre-processing step. However, additional information can be provided such as
a CNOlist (see makeCNOlist and readMIDAS) or information related to the 
pre-processing steps (compression, NONC nodes, expansion gates).
It can also be used to plot optimised model given the optimisation bitstring.}
\usage{
plotModel(model, cnolist=NULL, bString=NULL, indexIntegr=NA, signals=NULL, stimuli=NULL,
inhibitors=NULL, ncno=NULL, compressed=NULL, output="STDOUT",filename=NULL, 
graphvizParams=list())

}
\arguments{
  \item{model}{a model as returned by \code{readSif}. Alternatively, the filename can also be provided.}
  \item{cnolist}{output of makeCNOlist}
  \item{bString}{a sequence made of numbers between 0 and 1 of same length as
the one returned by the Genetic Algorithm (GA). This is a generalisation of
the bitString returned by the GA function: several bit strings can be
averaged and used}
  \item{indexIntegr}{additional indices to highlight some edge (optional).}
  \item{signals}{a list of nodes belonging to the signals class}
  \item{stimuli}{a list of nodes belonging to the stimuli class}
  \item{inhibitors}{a list of nodes belonging to the inhibitors class}
  \item{ncno}{a list of NCNO nodes.}
  \item{compressed}{a list of compressed nodes}
  \item{filename}{the filename (without extension) used to write the dot file}
  \item{output}{the type of output (PNG, PDF, SVG accepted)}
  \item{graphvizParams}{a list of optional arguments dedicated to Rgraphviz.
Accepted arguments are: \code{arrowsize} (default is 2); 
\code{size} (a string for the dot output default is "15,15");\code{fontsize}
(default is 22);\code{edgecolor} (default is "black").
}}




\details{
This function plots the model and also saves it in a dot file that can be processed
later on. However, you can also save the plot in PNG or PDF or SVG format (one at
a time).

The cnolist argument contains the signals/stimuli/inhibitors so if you provide
a cnolist there is no need to use these arguments. If you decide to use them
they will overwrite the contents of the cnolist argument.

optimRes is the output of gaBinary. One of its field is called bString and
contains a list of 0 and 1 (the optimisation is perfomed with a binary
procedure). This list of 0 and 1 is then used to plot or not the edges of the
model. However, you can provide a bitString made of floats (e.g., average of 
several bitStrings). In such case, edges will appear in gray light or dark
according to the bistring value (between 0 and 1).

}
\value{
a graph representation of the model
\item{graph$g}{A graph representation of the model}
\item{graph$attrs}{graph attributes}
\item{graph$nodeAttrs}{nodes attributes}
\item{graph$edgeAttrs}{edges attributes}
\item{graph$clusters}{clusters of nodes}


}
%\references{
%J. Saez-Rodriguez, L. G. Alexopoulos, J. Epperlein, R. Samaga, D. A. Lauffenburger, S. Klamt and P. K. Sorger. Discrete logic modeling as a means to link protein signaling networks with functional analysis of
%mammalian signal transduction, Molecular Systems Biology, 5:331, 2009.
%}
\author{
T. Cokelaer
}
\seealso{
\link{readMIDAS}, \link{readSif}, \link{makeCNOlist}, \link{writeNetwork},
\link{writeDot}, and \link{gaBinaryT1}
}

\note{This function depends on the Rgraphviz package.}

\examples{
    \dontrun{library(CellNOptR)
    data(CNOlistToy,package="CellNOptR")
    data(ToyModel,package="CellNOptR")
    res<-plotModel(ToyModel, cnolist=CNOlistToy, compressed=c("TRAF6", "p38"))
    }
}


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
plotModel(model, CNOlist=NULL, bString=NULL, indexIntegr=NULL, signals=NULL, 
stimuli=NULL, inhibitors=NULL, NCNO=NULL, compressed=NULL, output="STDOUT",
filename=NULL, graphvizParams=list(), show=TRUE, remove_dot=TRUE)
}

\arguments{
    \item{model}{
      a model as returned by \code{readSIF}. Alternatively, the filename can
      also be provided.
    }
    \item{CNOlist}{output of makeCNOlist}
    \item{bString}{
        a sequence made of numbers between 0 and 1 of same length as
        the one returned by the Genetic Algorithm (GA). This is a generalisation
        of the bitString returned by the GA function: several bit strings can be
        averaged and used.
    }
    \item{indexIntegr}{
        additional indices to highlight some edge (optional).
    }
    \item{signals}{
        a list of nodes belonging to the signals class
    }
    \item{stimuli}{
        a list of nodes belonging to the stimuli class
    }
    \item{inhibitors}{
        a list of nodes belonging to the inhibitors class
    }
    \item{NCNO}{
        a list of NCNO nodes.
    }
    \item{compressed}{
        a list of compressed nodes
    }
    \item{filename}{
        the filename (without extension) used to write the dot file
    }
    \item{output}{
        the type of output (PNG, PDF, SVG accepted)
    }
    \item{graphvizParams}{a list of optional arguments dedicated to Rgraphviz to
        tune the layout:
        \itemize{
        \item{arrowsize}{ default is 2 }
        \item{size}{a string for the size of the dot output; default is "15,15"}
        \item{fontsize}{ default is 22}
        \item{edgecolor}{ default is "black"}
        \item{nodeLabels}{ overwrit node label with a list of proper length.}
        \item{nodeWidth}{ default is 2}
        \item{nodeHeight}{ default is 1}
        \item{viewEmptyEdges}{  default is TRUE}
        \item{mode}{  can be 'classic' or 'sbgn' (default). The difference appears in the and gate.}
        \item{andWidth}{ = 0.2 in 'classic' mode and 0.5 in 'sbgn' mode}
        \item{andHeight}{ = 0.2 in 'classic' mode and 0.5 in 'sbgn' mode}
        }
    }

  \item{show}{ show the plot (default is True)}

  \item{remove_dot}{remove the dot file that has been created .default if
False)}
}




\details{
This function plots the model and also saves it in a dot file that can be processed
later on. However, you can also save the plot in PNG or PDF or SVG format (one at
a time).

The CNOlist argument contains the signals/stimuli/inhibitors so if you provide
a CNOlist there is no need to use these arguments. If you decide to use them
they will overwrite the contents of the CNOlist argument.

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
\author{
T. Cokelaer
}
\seealso{
\link{readMIDAS}, \link{readSIF}, \link{makeCNOlist}, \link{writeNetwork},
\link{writeDot}, \link{gaBinaryT1}
}

\note{This function depends on the Rgraphviz package.}

\examples{
    data(CNOlistToy,package="CellNOptR")
    data(ToyModel,package="CellNOptR")
    res<-plotModel(ToyModel, CNOlist=CNOlistToy, compressed=c("TRAF6", "p38"), 
        graphvizParams=list(mode="classic", fontsize=30))
    
}


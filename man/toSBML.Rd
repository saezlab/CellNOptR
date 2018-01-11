\name{toSBML}
\alias{toSBML}
\title{
Export the network to SBML-qual format
}
\description{
Export a Boolean network <network> to an sbml-qual file <fileName>.
This file can then be read in using other software that supports SBMLqual standards.

The function also takes a bit string as input. 
It cuts the model according to the values in bitstrings and write the new model object to SBMLqual.

}
\usage{
toSBML(network, file, bitString = c(rep(1,length(network$reacID))),version=c("standard","cellnopt"))
}
\arguments{
	\item{network}{a valid CellNOptR network model created by e.g. readSIF()}
  \item{file}{a valid filename to save the SBMLqual model}
  \item{bitString}{optional vector of binary values (for example the resulted bitString from optimisation) to cut the unnecessary interactions from the network before exporting.}
   \item{version}{ defines the format of SBMLqual file, read details.}
}
\details{
	\code{version = "standard"} exports one transition block for each node of the network. This format is the SBMLqual standard, that can be imported then with other softwares
	
	\code{version = "cellnopt"} the exported file follows a simplified syntax, where each edge of the network is transformed to a transition block in the SBMLqual file. Can be later imported to CellNOptR again.
}

\author{
Francesco Ceccarelli
}



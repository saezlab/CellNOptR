\name{model2sif}
\alias{model2sif}

\title{
  Convert a model object in sif format
}

\description{
  This function receives as input a model object and converts it to the cytoscape sif format. It can be used to convert either the whole model (before or after pre-processing) or the optimized one (if the corresponding bitString in provided). The resulting table can also be saved in a sif file.
}

\usage{
model2sif(model,optimRes=NA,writeSif=FALSE, filename="Model")

}


\arguments{
  \item{model}{
the model as generated using \link{readSIF}
}
  \item{optimRes}{
the output of the optimisation (as obtained using \link{gaBinaryT1}), default set to NA the whole model in converted 
}
  \item{writeSif}{
it writeSif=FALSE (default) the sif file is not saved. If writeSif=TRUE it is saved.
}
  \item{filename}{
the name of the sif file saved if writeSif=TRUE. Default is Model.sif.
}

}


\details{
All links in the model are converted in sif format that is sourceNode-tab-sign-tab-targetNode. If there are ANDs, they are converted using dummy nodes called and# (e.g. A+B=C will be A-tab-1-tab-and1; B-tab-1-tab-and1; and1-tab-1-tab-C).
}

\value{
\item{sifFile}{a table with all the links in the model in the format sourceNode-tab-sign-tab-targetNode}
}

\author{
    F.Eduati
}

\seealso{
\link{graph2sif}, \link{sif2graph}, \link{readSIF}, 
}



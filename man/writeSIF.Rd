\name{writeSIF}
\alias{writeSIF}
\alias{writeSif}

\title{
Convert a model into a SIF format and save the result in a file.
}
\description{
This function takes as input a model (as created by e.g., read from a SIF data set with
readSIF function) and save it into a file. 
}
\usage{
writeSIF(model, filename,  overwrite = FALSE)
}
\arguments{
  \item{model}{the model
  }
  \item{filename}{the filename
  }
  \item{overwrite}{by default, do not overwrite a file.
  }


}

\author{
T Cokelaer
}
\examples{
tmpdir<-tempdir()
setwd(tmpdir)
cpfile<-dir(system.file("ToyModel",package="CellNOptR"),full=TRUE)
file.copy(from=cpfile,to=getwd(),overwrite=TRUE)
ToyModel<-readSIF(sifFile="ToyPKNMMB.sif")
writeSIF(ToyModel, "ToyPKNMMB_copy.sif")
}

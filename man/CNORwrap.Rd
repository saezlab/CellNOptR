\name{CNORwrap}
\alias{CNORwrap}
\title{
CNOR analysis wrapper
}
\description{
This function is a wrapper around the whole CNOR analysis, it performs the following steps:

 \enumerate{
  \item Plot the CNOlist;
  \item Checks data to model compatibility;
  \item Call the preprocessing function (Cut the nonc off the model, compress the model and expand the gates);
  \item Compute the residual error;
  \item Prepare for simulation;
  \item Optimisation T1 and T2 (optional);
  \item Plot simulated and experimental results;
  \item Plot the evolution of fit;
  \item Write the scaffold and PKN;
  \item Write the report
 }
}
\usage{
CNORwrap(paramsList=NA, data=NA, model=NA, name, namesData=NA, time=1,
compression=TRUE, expansion=TRUE, cutNONC=TRUE)
}

\arguments{

  \item{paramsList}{
    A list of parameters related to the Genetic Algorithm parameters:
   \enumerate{
     \item sizeFac: default to 1e-04;
     \item NAFac: default to 1;popSize: default to 50;
     \item pMutation: default to 0.5;
     \item maxTime: default to 60;
     \item maxGens: default to 500;
     \item stallGenMax: default to 100;
     \item selPress: default to 1.2;
     \item elitism: default to 5;
     \item relTol: default to 0.1;
     \item verbose: default to FALSE (default to true in the functions used by CNORwrap but CNORwrap sets them to false by default).
   }
  }
  and the Data and Model structure.

  If paramsList is not provided (NA), it is filled internally with the \link{defaultParameters} function.

  If Data and Model are not provided in \code{paramList}, the function looks for Data and Model arguments.

  If Data and Model are provided, the function overwrites the field \code{data} and \code{model} in \code{paramsList}.


  \item{data}{
a CNOlist structure, as created by \link{makeCNOlist}
}
  \item{model}{
a model structure, as created by \link{readSIF}.
}
 \item{name}{
a string that will be used to name the project and all graphs produced
}
 \item{namesData}{
a list with two elements:CNOlist and Model, each containing a string that is a
reference for the user to know which model/data set was used (it will be
included in the report). If not provided, the list is built automatically using
the name arguments.
}
 \item{time}{
either 1 or 2: Do you want to perform a one time point steady state optimisation
or a 2 time points pseudo steady state optimisation. By default this is set to 1.
}
\item{compression}{
compress the model (default TRUE)
}
\item{expansion}{
expand the gates (default TRUE)

}
\item{cutNONC}{
cut the NONC nodes off the model by (default TRUE)

}

}
\details{
If you do not provide a parameters list, you can provide only essential elements,
and all other parameters will be set to their default values. In this case, you
should set paramsList=NA, and provide the following fields: data, model, name, time.
}
\value{
This function does not return anything, it does the analysis, produces all the
plots and puts them in a folder that is in your working directory, and is called "Name".
}
\author{
C. Terfve
}
\examples{
#version with paramslist

tmpdir<-tempdir()
setwd(tmpdir)

data(CNOlistToy,package="CellNOptR")
data(ToyModel,package="CellNOptR")

pList = defaultParameters(CNOlistToy, ToyModel)
pList$maxGens = 5
pList$popSize = 5

CNORwrap(paramsList=pList, name="Toy",
    namesData=list(CNOlist="ToyData", model="ToyModel"))

\dontrun{
#version with default parameters

data(CNOlistToy,package="CellNOptR")
data(ToyModel,package="CellNOptR")

CNORwrap(name="Toy",
    namesData=list(CNOlist="ToyData",model="ToyModel"),
	data=CNOlistToy, model=ToyModel)
}
}

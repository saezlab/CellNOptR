\name{CNORwrap}
\alias{CNORwrap}
\title{
CNOR analysis wrapper
}
\description{
This function is a wrapper around the whole CNOR analysis, it performs the following steps: 1. Plot the CNOlist; 2. Checks data to model compatibility; 3. Find the indices, in the model, of the species that are inhibited/stimulated/measured (signals); 4. Find the indices of the non-osbervables/non-controllables (nonce); 5. Cut the nonc off the model; 6. Recompute the indices; 7. Compress the model;8. Recompute the indices; 9. Expand the gates; 10. Compute the residual error; 11. Prepare for simulation; 12. Optimisation; 13. Plot simulated and experimental results; 14. Plot the evolution of fit; 15. Write the scaffold and PKN; 16. Write the report
}
\usage{
CNORwrap(paramsList,Data,Model,Name,NamesData,Time=1)
}
\arguments{
  \item{paramsList}{
paramsList has entries: Data:the CNOlist, Model:the model; Parameters for the optimisation:sizeFac: default to 1e-04; NAFac: default to 1;PopSize: default to 50; Pmutation: default to 0.5; MaxTime: default to 60; maxGens: default to 500; StallGenMax: default to 100; SelPress: default to 1.2; elitism: default to 5; RelTol: default to 0.1; verbose: default to FALSE (default to true in the functions  used by CNORwrap but CNORwrap sets them to false by default).
}
  \item{Data}{
the CNOlist that contains the data that you will use
}
  \item{Model}{
the model that you want to optimise
}
 \item{Name}{
a string that will be used to name the project and all graphs produced
}
 \item{NamesData}{
a list with two elements:CNOlist and Model, each containing a string that is a reference for the user to know which model/data set was used (it will be included in the report)
}
 \item{Time}{
either 1 or 2: Do you want to perform a one time point steady state optimisation or a 2 time points pseudo steady state optimisation. By default this is set to 1.
}
}
\details{
If you do not provide a parameters list, you can provide only essential elements, and all other parameters will be set to their default values. In this case, you should set paramsList=NA, and provide the following fields: Data, Model, Name, Time.
}
\value{
This function does not return anything, it does the analysis, produces all the plots and puts them in a folder that is in your working directory, and is called "Name".
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

pList<-list(
	Data=CNOlistToy,
	Model=ToyModel,
	sizeFac = 1e-04, 
	NAFac = 1, 
	PopSize = 5, 
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

\dontrun{
#version with default parameters

dir.create("CNOR_analysis")
setwd("CNOR_analysis")

data(CNOlistToy,package="CellNOptR")
data(ToyModel,package="CellNOptR")

CNORwrap(
	paramsList=NA,
	Name="Toy",
	NamesData=list(CNOlist="ToyData",Model="ToyModel"),
	Data=CNOlistToy,
	Model=ToyModel)
}
}

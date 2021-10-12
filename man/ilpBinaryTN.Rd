\name{ilpBinaryTN}
\alias{ilpBinaryTN}
\title{
ILP method used to optimise a model
}
\description{
This function is the ilp method to be used to optimise a model by fitting to data for time point 2, that should follow optimisation based on time point 1.
}
\usage{
ilpBinaryTN(cnolist, 
                        model, 
                        sizeFac = 0.0001, 
                        mipGap = 0, 
                        relGap = 0, 
                        timelimit = 3600, 
                        cplexPath, 
                        method = "quadratic",
                        numSolutions = 100, 
                        limitPop = 500, 
                        poolIntensity = 0, 
                        poolReplace = 2,
                        timeIndices = c(1, 2))
}
\arguments{
  \item{cnolist}{
a CNOlist on which the score is based (based on valueSignals[[2]], i.e. data at
time 1)
}
  \item{model}{
a model structure, as created by \code{readSIF}, normally pre-processed but that
is not a requirement of this function
}
  \item{sizeFac}{
the scaling factor for the size term in the objective function, default to 0.0001
}
  \item{mipGap}{
the absolute tolerance on the gap between the best integer objective and the objective of the best node remaining. When this difference falls below the value of this parameter, the linear integer optimization is stopped. Default set to 0
}
  \item{relGap}{
the relative tolerance on the objective value for the solutions in the solution pool. Solutions that are worse (either greater in the case of a minimization, or less in the case of a maximization) than the incumbent solution by this measure are not kept in the solution pool. Default set to 0
}
  \item{timelimit}{
the maximum optimisation time in seconds, default set to 3600
}
  \item{cplexPath}{
the path where the cplex solver is stored. Default set to "~/Documents/cplex"
}
  \item{method}{
the method of writing the objective function (quadratic/linear). Default set to "quadratic"
}
  \item{numSolutions}{
the number of solutions to save
}
  \item{limitPop}{
the number of solutions to be generated. Default set to 500
}
  \item{poolIntensity}{
the Intensity of solution searching. Default set to 4
}
  \item{timeIndices}{
the time indeces to optimize. Default set to timeIndices=c(1, 2)
}

  \item{poolReplace}{
 pool replacement strategy, consult CPLEX manual for details. 
}

}
\value{
This function returns a list with elements:
\item{bitstringILPAll}{the list of all optimal bitstrings identified}
\item{bScore}{the best score for each set of bitstrings}
\item{time_cplex_only}{the time it took for cplex to solve the problem}
\item{total_time}{the total time for the pipeline to run (writing problem + solving problem + retrieving solutions)}
\item{stringsTolScores}{the scores of the above-mentioned strings}
}
\author{
E Gjerga, H Koch
}
\references{
Alexander Mitsos, Ioannis N. Melas, Paraskeuas Siminelakis, Aikaterini D. Chairakaki, Julio Saez-Rodriguez, and Leonidas G. Alexopoulos. Identifying Drug Effects via Pathway Alterations using an Integer Linear Programming Optimization Formulation on Phosphoproteomic Data. PLoS Comput Biol. 2009 Dec; 5(12): e1000591.
}
\examples{


# Toy Exampple
data("ToyModel", package="CellNOptR")
data("CNOlistToy", package="CellNOptR")
pknmodel = ToyModel
cnolist = CNOlist(CNOlistToy)
model = preprocessing(data = cnolist, model = pknmodel, compression = TRUE, expansion = TRUE)
plotModel(model = model, CNOlist = cnolist)

# Training to data - ILP
\dontrun{
resILP = ilpBinaryTN(cnolist = cnolist, model = model)
}
}

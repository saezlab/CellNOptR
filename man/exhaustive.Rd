\name{exhaustive}
\alias{exhaustive}
\title{Exhaustive search over the optimisation of a PKN model on MIDAS data.}

\description{
This function performs an exhaustive search of the parameter space tring all the
solutions. It is used internally by the genetic algorithm when a small model has
to be optimised and the number of solutions to try is smaller than the number of
iterations that the Genetic Algorithm will perform.
}
\usage{
exhaustive(CNOlist, model, shuffle=FALSE, Nmax=NULL, verbose=TRUE, sizeFac =
0.0001, NAFac = 1, relTol=0.1, timeIndex=2)
}
\arguments{
  \item{CNOlist}{
a CNOlist on which the score is based (based on valueSignals[[2]], i.e. data at
time 1)
}
  \item{model}{
a model structure, as created by \code{readSIF}, normally pre-processed but that
is not a requirement of this function
}
  \item{shuffle}{
The list of bitstrings is set up arbitrarely. You may want to shuffle it.
}
  \item{Nmax}{
The total number of computation will be 2 to the power N, where N is the size of
the model (ReacID field). The total number of computation can be large. You may
want to set a maximumn number of computation using Nmax. 
}
  \item{sizeFac}{
the scaling factor for the size term in the objective function, default to 0.0001
}
  \item{NAFac}{
the scaling factor for the NA term in the objective function, default to 1
}
  \item{relTol}{
the relative tolerance for the best bitstring reported by the genetic algorithm,
i.e., how different from the best solution, default set to 0.1
Not yet implemented.
}
  \item{verbose}{
logical (default to TRUE) do you want the statistics of each generation to be printed on the screen?
}

  \item{timeIndex}{the index of the time point to optimize. Must be greater or
equal to 2 (1 corresponds to time=0). Must be less than the number of time
points. Default is 2. 
    }


}
\value{
This function returns a list with elements:
\item{bString}{the best bitstring}
\item{bScore}{the best score}
\item{all_scores}{all scores that have been computed}
\item{results}{ a matrix with columns "Generation","Best_score","Best_bitString","Stall_Generation","Avg_Score_Gen","Best_score_Gen","Best_bit_Gen","Iter_time"}
\item{stringsTol}{the bitstrings whose scores are within the tolerance}
\item{stringsTolScores}{the scores of the above-mentioned strings}

Note that the field results, is not yet populated but maybe in the future.
}
\author{
T. Cokelaer
}
\seealso{
\link{gaBinaryT1}
}
\examples{
data(CNOlistToy,package="CellNOptR")
data(ToyModel,package="CellNOptR")

#pre-process model

model = preprocessing(CNOlistToy, ToyModel)

#optimise

results <-exhaustive(
	CNOlist=CNOlistToy,
	model=model,
    shuffle=TRUE,
    Nmax=1000,
	verbose=FALSE)
}

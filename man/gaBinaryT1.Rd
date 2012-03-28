\name{gaBinaryT1}
\alias{gaBinaryT1}
\title{
Genetic algorithm used to optimise a model
}
\description{
This function is the genetic algorithm to be used to optimise a model by fitting to data containing one time point.
}
\usage{
gaBinaryT1(CNOlist, Model, SimList, indexList, sizeFac = 1e-04, NAFac = 1, initBstring, PopSize = 50, Pmutation = 0.5, MaxTime = 60, maxGens = 500, StallGenMax = 100, SelPress = 1.2, elitism = 5, RelTol = 0.1, verbose=TRUE)
}
\arguments{
  \item{CNOlist}{
a CNOlist on which the score is based (based on valueSignals[[2]], i.e. data at t1)
}
  \item{Model}{
a Model list
}
  \item{SimList}{
a list that contains additional fields for the simulator, as created by prep4Sim applied to the model above
}
  \item{indexList}{
a list of indexes of species stimulated/inhibited/signals, as produced by indexfinder applied on the model and CNOlist above
}
  \item{sizeFac}{
the scaling factor for the size term in the objective function, default to 0.0001
}
  \item{NAFac}{
the scaling factor for the NA term in the objective function, default to 1
}
  \item{initBstring}{
an initial bitsring to be tested, should be of the same size as the number of reactions in the model above
}
  \item{PopSize}{
the population size for the genetic algorithm, default set to 50
}
  \item{Pmutation}{
the mutation probability for the genetic algorithm, default set to 0.5
}
  \item{MaxTime}{
the maximum optimisation time in seconds, default set to 60
}
  \item{maxGens}{
the maximum number of generations in the genetic algorithm, default set to 500
}
  \item{StallGenMax}{
the maximum number of stall generations in the genetic algorithm, default to 100
}
  \item{SelPress}{
the selective pressure in the genetic algorithm, default set to 1.2
}
  \item{elitism}{
the number of best individuals that are propagated to the next generation in the genetic algorithm, default set to 5
}
  \item{RelTol}{
the relative tolerance for the best bitstring reported by the genetic algorithm,
i.e., how different from the best solution, default set to 0.1
}
  \item{verbose}{
logical (default to TRUE) do you want the statistics of each generation to be printed on the screen?
}
}
\details{
The whole procedure is described in details in Saez-Rodriguez et al. (2009). The basic principle is that at each generation, the algorithm evaluates a population of models based on excluding or including some gates in the initial pre-processed model (this is encoded in a bitstring with contains 0/1 entries for each gate).  The population is then evolved based on the results of the evaluation of these networks, where the evaluation is obtained by simulating the model (to steady state) under the various conditions present in the data, and then computing the squared deviation from the data, to which a penalty is added for size of the model and for species in the model that do not reach steady state.
}
\value{
This function returns a list with elements:
\item{bString}{the best bitstring}
\item{Results}{ a matrix with columns "Generation","Best_score","Best_bitString","Stall_Generation","Avg_Score_Gen","Best_score_Gen","Best_bit_Gen","Iter_time"}
\item{StringsTol}{the bitstrings whose scores are within the tolerance}
\item{StringsTolScores}{the scores of the above-mentioned strings}
}
\references{
J. Saez-Rodriguez, L. G. Alexopoulos, J. Epperlein, R. Samaga, D. A. Lauffenburger, S. Klamt and P. K. Sorger. Discrete logic modeling as a means to link protein signaling networks with functional analysis of
mammalian signal transduction, Molecular Systems Biology, 5:331, 2009.
}
\author{
C. Terfve
}
\seealso{
GetFit, prep4Sim, indexFinder, simulatorT1
}
\examples{
data(CNOlistToy,package="CellNOptR")
data(ToyModel,package="CellNOptR")

#pre-process model

indicesToy<-indexFinder(CNOlistToy,ToyModel,verbose=FALSE)
ToyNCNOindices<-findNONC(ToyModel,indicesToy,verbose=FALSE)
ToyNCNOcut<-cutNONC(ToyModel,ToyNCNOindices)
indicesToyNCNOcut<-indexFinder(CNOlistToy,ToyNCNOcut)
ToyNCNOcutComp<-compressModel(ToyNCNOcut,indicesToyNCNOcut)
indicesToyNCNOcutComp<-indexFinder(CNOlistToy,ToyNCNOcutComp)
ToyNCNOcutCompExp<-expandGates(ToyNCNOcutComp)

#optimise

ToyFields4Sim<-prep4Sim(ToyNCNOcutCompExp)
initBstring<-rep(1,length(ToyNCNOcutCompExp$reacID))
ToyT1opt<-gaBinaryT1(
	CNOlist=CNOlistToy,
	Model=ToyNCNOcutCompExp,
	SimList=ToyFields4Sim,
	indexList=indicesToyNCNOcutComp,
	initBstring=initBstring,
	maxGens=2,
	PopSize=5,
	verbose=FALSE)
}

\name{gaBinaryT2}
\alias{gaBinaryT2}
\title{
Genetic algorithm for time point 2
}
\description{
This is the genetic algorithm for time point 2, that should follow optimisation based on time point 1.
}
\usage{
gaBinaryT2(CNOlist, Model, SimList, indexList, bStringT1, SimResT1, sizeFac = 1e-04, NAFac = 1, PopSize = 50, Pmutation = 0.5, MaxTime = 60, maxGens = 500, StallGenMax = 100, SelPress = 1.2, elitism = 5, RelTol = 0.1, verbose=TRUE)
}
\arguments{
  \item{CNOlist}{
a CNOlist on which the score is based (based on valueSignals[[3]], i.e. data at  t2)
}
  \item{Model}{
a Model list
}
  \item{SimList}{
a list that contains additional fields for the simulator, as created by prep4Sim applied to the model above
}
  \item{indexList}{
a list of indexes of species stimulated/inhibited/signals, as produced by indexfinder applied on teh model and CNOlist above
}
  \item{bStringT1}{
the optimal bitstring from optimisation at time 1 (i.e. a vector of 0s and 1s)
}
  \item{SimResT1}{
the results of optimisation at t1 (i.e. the whole list as output by gabinaryT1)
}
  \item{sizeFac}{
the scaling factor for the size term in the objective function, default to 0.0001
}
  \item{NAFac}{
the scaling factor for the NA term in the objective function, default to 1
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
the number of best individuals that are propagated to the next generation in the gen. al. default set to 5
}
  \item{RelTol}{
the relative tolerance for the best bitstring reported by the genetic algorithm, i.e.how different from the best solution can solutions be to be reported as well, default set to 0.1
}
  \item{verbose}{
logical (default to TRUE) do you want the statistics of each generation to be printed on the screen
}
}
\details{
This function takes in the same input as the T1 ga, but in addition it takes in the bitstring optimised for T1, and does not take an initial bitstring. It also takes a matrix SimResT1, which is the results of simulatorT1 ran on the model optimised for T1. Be aware that the bitString that this function returns is one that only includes the bits that it actually looks at, i.e. the bits that were 0 in the bStringT1
}
\value{
This function returns a list with elements:
\item{bString}{the best bitstring}
\item{Results}{ a matrix with columns "Generation","Best_score","Best_bitString","Stall_Generation","Avg_Score_Gen","Best_score_Gen","Best_bit_Gen","Iter_time"}
\item{StringsTol}{the bitstrings whose scores are within the tolerance}
\item{StringsTolScores}{the scores of the above-mentionned strings}
}
\author{
C.Terfve
}
\seealso{
GetFit, prep4Sim, indexFinder, simulatorT1, simulatorT2, gaBinaryT2
}
\examples{
#load data

data(CNOlistToy2,package="CellNOptR")
data(ToyModel2,package="CellNOptR")

#pre-process model

checkSignals(CNOlistToy2,ToyModel2)
indexesToy2<-indexFinder(CNOlistToy2,ToyModel2,verbose=FALSE)
ToyNCNOindexes2<-findNONC(ToyModel2,indexesToy2,verbose=FALSE)
ToyNCNOcut2<-cutNONC(ToyModel2,ToyNCNOindexes2)
indexesToyNCNOcut2<-indexFinder(CNOlistToy2,ToyNCNOcut2)
ToyNCNOcutComp2<-compressModel(ToyNCNOcut2,indexesToyNCNOcut2)
indexesToyNCNOcutComp2<-indexFinder(CNOlistToy2,ToyNCNOcutComp2)
ToyNCNOcutCompExp2<-expandGates(ToyNCNOcutComp2)

#optimise t1
ToyFields4Sim2<-prep4Sim(ToyNCNOcutCompExp2)
initBstring2<-rep(1,length(ToyNCNOcutCompExp2$reacID))
ToyT1opt2<-gaBinaryT1(
	CNOlist=CNOlistToy2,
	Model=ToyNCNOcutCompExp2,
	SimList=ToyFields4Sim2,
	indexList=indexesToyNCNOcutComp2,
	initBstring=initBstring2,
	maxGens=2, 
	PopSize = 5,
	verbose=FALSE)

#Optimise T2
SimToyT12<-simulateT1(
	CNOlist=CNOlistToy2,
	Model=ToyNCNOcutCompExp2,
	bStringT1=ToyT1opt2$bString,
	SimList=ToyFields4Sim2,
	indexList=indexesToyNCNOcutComp2)
ToyT1opt2T2<-gaBinaryT2(
	CNOlist=CNOlistToy2,
	Model=ToyNCNOcutCompExp2,
	SimList=ToyFields4Sim2,
	indexList=indexesToyNCNOcutComp2,
	bStringT1=ToyT1opt2$bString,
	SimResT1=SimToyT12,
	maxGens=2, 
	PopSize = 5,
	verbose=FALSE)
}

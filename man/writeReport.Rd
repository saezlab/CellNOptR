\name{writeReport}
\alias{writeReport}
\title{
Write a report of a CellNOptR analysis
}
\description{
This function writes a short report of a CellNOptR analysis in an html page, that is linked to the various graphs produced
}
\usage{
writeReport(modelOriginal, modelOpt, optimResT1, optimResT2, CNOlist, directory,
    namesFiles = list(dataPlot = NA, evolFitT1=NA, evolFitT2=NA, simResultsT1=NA,
    simResultsT2=NA, scaffold=NA, scaffoldDot=NA, tscaffold=NA, wscaffold=NA,
    PKN=NA, PKNdot=NA, wPKN=NA,nPKN=NA), namesData = list(CNOlist=NA, model=NA),
    resE=NULL)
}


\arguments{
  \item{modelOriginal}{
the original previous knowledge network (i.e. model that you loaded) in a model list format
}
  \item{modelOpt}{
the model that was actually used for optimisation (i.e. the scaffold network, after compression and expansion) in a model list format
}
  \item{optimResT1}{
the results of the optimisation at t1, as output by gabinaryT1
}
  \item{optimResT2}{
the results of the optimisation at t2, as output by gabinaryT2. Set this to NA if you have performed a one time point optimisation.
}
  \item{CNOlist}{
a CNOlist
}
  \item{directory}{
the name of a new directory that will be created, where your results will be moved
}
  \item{namesFiles}{
 a list of the names of the files that should have been created. Depending on
whether a t2 optimisation was performed or not, all or some of the following
fields are expected: dataPlot, evolFitT1, evolFitT2, simResultsT1, simResultsT2,
scaffoldDot, scaffold, tscaffold, wscaffold, PKN, PKNdot, wPKN, nPKN.
}
  \item{namesData}{
a list with fields \code{$CNOlist} and \code{$model} that contain strings that
are meaningful identifiers of your data and previous knowledge network (for your
own record).
}
  \item{resE}{
a vector with named entries t1, t2 t1andt2, as produced by the function
\code{ResidualError}, that contains the residual error associated with the
discretisation of the data. Since version 1.3.29, there is no need to provide
this argument. Kept for back compatibility.
}
}
\details{
Future versions of this function might directly write and compile a tex file.
}
\value{
This function produces a directory and moves all the files of namesFiles to it,
then it creates an html report that contains infos about the optimisation
process.
}
\author{
C. Terfve
}
\seealso{
\link{writeNetwork}, \link{writeScaffold}
}
\examples{
tmpdir<-tempdir()
setwd(tmpdir)

#load data

data(CNOlistToy,package="CellNOptR")
data(ToyModel,package="CellNOptR")

#pre-process model (partial)

indicesToy<-indexFinder(CNOlistToy,ToyModel,verbose=TRUE)
ToyNCNOcutComp<-compressModel(ToyModel,indicesToy)
indicesToyNCNOcutComp<-indexFinder(CNOlistToy,ToyNCNOcutComp)
ToyNCNOcutCompExp<-expandGates(ToyNCNOcutComp)

#optimise

ToyFields4Sim<-prep4sim(ToyNCNOcutCompExp)
initBstring<-rep(1,length(ToyNCNOcutCompExp$reacID))
ToyT1opt<-gaBinaryT1(
	CNOlist=CNOlistToy,
	model=ToyNCNOcutCompExp,
	initBstring=initBstring, 
	maxGens=2,
	popSize=5,
	verbose=TRUE)

#write report

namesFilesToy<-list(
	dataPlot=NA,
	evolFitT1=NA,
	evolFitT2=NA,
	simResultsT1=NA,
	simResultsT2=NA,
	scaffold=NA,
	scaffoldDot=NA,
	tscaffold=NA,
	wscaffold=NA,
	PKN=NA,
	PKNdot=NA,
	wPKN=NA,
	nPKN=NA)

writeReport(
	modelOriginal=ToyModel,
	modelOpt=ToyNCNOcutCompExp,
	optimResT1=ToyT1opt,
	optimResT2=NA,
	CNOlist=CNOlistToy,
	directory="testToy",
	namesFiles=namesFilesToy,
	namesData=list(CNOlist="Toy",model="ToyModel"),
	resE=NA)
}

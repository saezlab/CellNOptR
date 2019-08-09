\name{invokeCPLEX}
\alias{invokeCPLEX}
\alias{invokeCPLEX}

\title{
Solving the ILP problem with CPLEX.
}
\description{
This function takes as an input the file name where we write the ILP formulation together with CPLEX parameters
and then solves the ILP problem. 
}
\usage{
invokeCPLEX <- function(inputFileName, 
                        outputFileName, 
                        mipGap=mipGap, 
                        relGap = relGap, 
                        timelimit=timelimit, 
                        cplexPath = cplexPath, 
                        numSolutions = numSolutions, 
                        limitPop = limitPop, 
                        poolIntensity = poolIntensity, 
                        poolReplace = poolReplace)
}
\arguments{
  \item{inputFileName}{the file name where the cplex ilp problem is stored
  }
  \item{outputFileName}{the file name where to store the cplex result
  }
  \item{mipGap}{the mipgap
  }
  \item{relGap}{the relGap
  }
  \item{timelimit}{the timelimit
  }
  \item{cplexPath}{the cplex solver path
  }
  \item{numSolutions}{the number of desired solutions
  }
  \item{limitPop}{the limitPop
  }
  \item{poolIntensity}{the poolIntensity
  }
  \item{poolReplace}{the poolReplace
  }


}

\author{
E Gjerga, H Koch
}
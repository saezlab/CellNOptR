\name{writeFile}
\alias{writeFile}
\alias{writeFile}

\title{
Writing the ILP problem.
}
\description{
This function takes as input the objective function. constraints, bounds and the solver path in order to generate a file containing the ILP problem. 
}
\usage{
writeFile <- function(objectiveFunction, constraints, bounds, binaries,
                      cplexPath)
}
\arguments{
  \item{objectiveFunction}{the objective function of the ILP problem
  }
  \item{constraints}{the set of constraints of the ILP problem
  }
  \item{bounds}{the set of bounds for each integer variable
  }
  \item{binaries}{the set of binary variables
  }
  \item{cplexPath}{the path of the cplex solver
  }


}
\value{
Writing the LP file which contains the ILP problem for CellNOpt.
}

\author{
E Gjerga, H Koch
}
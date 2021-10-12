\name{writeFile}
\alias{writeFile}

\title{
Writing the ILP problem.
}
\description{
This function takes as input the objective function. constraints, bounds and the solver path in order to generate a file containing the ILP problem. 
}
\usage{
writeFile(objectiveFunction,
                      constraints,
                      bounds,
                      binaries)
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

}

\author{
E Gjerga, H Koch
}
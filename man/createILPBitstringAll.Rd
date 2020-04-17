\name{createILPBitstringAll}
\alias{createILPBitstringAll}
\alias{createILPBitstringAll}

\title{
Reading the optimal solutions as bitstrings.
}
\description{
This function takes as the cplex optimization results and the variables assigned to each interaction in the ILP formulation
in order to read the optimal bitstring. 
}
\usage{
createILPBitstringAll<- function(cplexSolutionFileName, y_vector,
                                 binary_variables)
}
\arguments{
  \item{cplexSolutionFileName}{the file name where the cplex results are stored
  }
  \item{y_vector}{the variables for each interaction in the PKN
  }
   \item{binary_variables}{the binary variables of the ILP formulation
  }


}
\value{
Creating the binary bitstring from the optimal solution from the ILP analysis.
}
\author{
E Gjerga, H Koch
}
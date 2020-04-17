\name{createAndRunILP}
\alias{createAndRunILP}
\alias{createAndRunILP}
\title{
Creating and running the ILP problem.
}
\description{
This function takes as an input the cno inputs (model + data) together with CPLEX parameters
and then solves the ILP problem. 
}
\usage{
createAndRunILP <- function(model = model, midas = midas, cnolist = cnolist, 
                            accountForModelSize = accountForModelSize, 
                            sizeFac = sizeFac, source_path = source_path, 
                            mipGap = mipGap, relGap = relGap, 
                            timelimit = timelimit, cplexPath = cplexPath, 
                            method = method, numSolutions = numSolutions, 
                            limitPop = limitPop, poolIntensity = poolIntensity, 
                            poolReplace = poolReplace)
}
\arguments{
  \item{model}{the model
  }
  \item{midas}{the midas table
  }
  \item{cnolist}{the cnolist object
  }
  \item{accountForModelSize}{the verbose parameter whether to account for model size
  }
  \item{sizeFac}{the size penalty factor
  }
  \item{source_path}{the source path
  }
  \item{mipGap}{the mipgap
  }
  \item{relGap}{the relGap
  }
  \item{timelimit}{the timelimit
  }
  \item{cplexPath}{the cplex solver path
  }
  \item{method}{the optimization method (quadratic/linear)
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
\value{
Creating and running the ILP problem for CellNOpt.
}
\author{
E Gjerga, H Koch
}
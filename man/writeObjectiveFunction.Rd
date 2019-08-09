\name{writeObjectiveFunction}
\alias{writeObjectiveFunction}
\alias{writeObjectiveFunction}

\title{
Writing the objective function for the ILP implementation of CellNOptR.
}
\description{
This function takes as input the integer variables assigned to each of the elements in the network and data in a CNOlist object
in order to produce the objective function needed to be optimized by the ILP solver. 
}
\usage{
writeObjectiveFunction <- function(model, 
                                   midasExperimentPart, 
                                   y_vector=y_vector,
                                   accountForModelSize = TRUE, 
                                   sizeFac = .000001, 
                                   meansOfMeasurements_at_t0, 
                                   method = "quadratic" )
}
\arguments{
  \item{model}{the model
  }
  \item{midasExperimentPart}{the experimental part of CNOlist where data is stored
  }
  \item{y_vector}{the variables for each interaction in the PKN
  }
  \item{accountForModelSize}{the verbose variable if we wish to apply the size penalty
  }
  \item{sizeFac}{the size penalty factor
  }
  \item{meansOfMeasurements_at_t0}{the means of measurements at time-point 0
  }
  \item{method}{the method of defining the objective function ("quadratic/linear")
  }


}

\author{
E Gjerga, H Koch
}
\references{
Alexander Mitsos, Ioannis N. Melas, Paraskeuas Siminelakis, Aikaterini D. Chairakaki, Julio Saez-Rodriguez, and Leonidas G. Alexopoulos. Identifying Drug Effects via Pathway Alterations using an Integer Linear Programming Optimization Formulation on Phosphoproteomic Data. PLoS Comput Biol. 2009 Dec; 5(12): e1000591.
}

\name{write_constraints}
\alias{write_constraints}
\alias{write_constraints}

\title{
Writing the set of constraints for the ILP implementation of CellNOptR.
}
\description{
This function takes as input the integer variables assigned to each of the elements in the network and data in a CNOlist object
in order to produce the set of constraints needed to be optimized by the ILP solver. 
}
\usage{
write_constraints(model,
                              midasExperimentPart,
                              midasTreatmentPart,
                              reaction_sets,
                              y_vector, 
                              midas, 
                              binary_variables)
}
\arguments{
  \item{model}{the model
  }
  \item{midasExperimentPart}{the experimental part of CNOlist where data is stored
  }
  \item{midasTreatmentPart}{the treatment part of CNOlist where data is stored
  }
  \item{reaction_sets}{the set of reactions
  }
  \item{y_vector}{the variables for each interaction in the PKN
  }
  \item{midas}{the midas table
  }
  \item{binary_variables}{the set of binary variables
  }


}

\author{
E Gjerga, H Koch
}
\references{
Alexander Mitsos, Ioannis N. Melas, Paraskeuas Siminelakis, Aikaterini D. Chairakaki, Julio Saez-Rodriguez, and Leonidas G. Alexopoulos. Identifying Drug Effects via Pathway Alterations using an Integer Linear Programming Optimization Formulation on Phosphoproteomic Data. PLoS Comput Biol. 2009 Dec; 5(12): e1000591.
}

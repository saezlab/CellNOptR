\name{write_bounds}
\alias{write_bounds}
\alias{write_bounds}

\title{
Writing the set of boundaries for each integer variable for the ILP implementation of CellNOptR.
}
\description{
This function takes as input the integer variables assigned to each of the elements in the network and data in a CNOlist object
in order to produce the set boundaries for each integer variable. 
}
\usage{
write_bounds <- function(model,
                         midasTreatmentPart,
                         y_vector,
                         binary_variables)
}
\arguments{
  \item{model}{the model
  }
  \item{midasTreatmentPart}{the treatment part of CNOlist where data is stored
  }
  \item{y_vector}{the variables for each interaction in the PKN
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

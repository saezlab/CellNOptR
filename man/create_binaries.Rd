\name{create_binaries}
\alias{create_binaries}
\alias{create_binaries}

\title{
Defining the set of binary variables for the ILP implementation of CellNOptR.
}
\description{
This function takes as input the model and data in a CNOlist object
in order to define the set of binary variablese. 
}
\usage{
create_binaries(model, 
                            midas, 
                            numberOfExperiments, 
                            y_vector)
}
\arguments{
  \item{model}{the model
  }
  \item{midas}{the midas table
  }
  \item{numberOfExperiments}{the number of experimental conditions
  }
  \item{y_vector}{the variables for each interaction in the PKN
  }


}

\author{
E Gjerga, H Koch
}
\references{
Alexander Mitsos, Ioannis N. Melas, Paraskeuas Siminelakis, Aikaterini D. Chairakaki, Julio Saez-Rodriguez, and Leonidas G. Alexopoulos. Identifying Drug Effects via Pathway Alterations using an Integer Linear Programming Optimization Formulation on Phosphoproteomic Data. PLoS Comput Biol. 2009 Dec; 5(12): e1000591.
}

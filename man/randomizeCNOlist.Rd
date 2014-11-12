\name{randomizeCNOlist}
\alias{randomizeCNOlist}
\title{
add noise to the data contained in a CNOlist.
}
\description{
This function takes in a CNOlist and does the normalisation of the data between
0 and 1, according to two different procedures (see details).
}
\usage{
randomizeCNOlist(cnolist, sd=0.1, minValue=0, maxValue=1, mode="gaussian")
}


\arguments{
  \item{cnolist}{
a CNOlist
}
  \item{sd}{standard deviation to be used when adding gaussian noise. Not used
if mode is uniform.
}
  \item{minValue}{When adding Gaussian noise, the result may be below the
minValue(default 0). If so, the value is set to minValue.
}
  \item{maxValue}{When adding Gaussian noise, the result may be above the
maxValue(default 1). If so, the value is set to maxValue.
}
  \item{mode}{The mode can be either 'gaussian', 'shuffle' or uniform'. In gaussian mode, a
gaussian noise is added to the data. The mean parameter is the
data and the standard deviation is defined by the sd parameter. In uniform mode,
the data is simply replaced by values taken from a uniform distribution between
0 and 1. In 'shuffle' mode all rows and columns are shuffled.
}

}

\value{
a noisy CNOlist
}
\author{
T. Cokelaer
}
\examples{
data(CNOlistToyMMB, package="CellNOptR")
cnolist = CNOlistToyMMB
cnolist2 = randomizeCNOlist(cnolist, mode="uniform")

# a method called randomize is available in the CNOlist class so you could type:
cnolist2 = randomize(cnolist, mode="uniform")
}

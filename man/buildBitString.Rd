\name{buildBitString}
\alias{buildBitString}
\title{
Build the final bit string vector and times vector based on a list of optimised bit strings at different
time points.
}
\description{
This function takes as input a list of vectors (can be only one). Each vector
represents an optimised bit string at a different time point (as returned by 
the gaBinary functions). The first optimised bit string have the same length as
the model to be optimised. The length of the following vectors corresponds to
the number of zeros found in the previous bitstring. For instance,
the following list of bit strings \code{bStrings = list(c(1,1,1,0,0,0),
c(1,0,0))} is correct whereas \code{bStrings = list(c(1,1,1,0,0,0), c(1,0))} is
incorrect.


This function is used internally by \link{computeScoreTN} and
\link{simulateTN}. It should not be used by a user in principle.
However, it may be useful for post processing.

New in version 1.3.28.
}
\usage{
buildBitString(bStrings)
}
\arguments{
    \item{bStrings}{
        a list of bit strings as returned by the optimisation at different time points.
    }
}


\value{
This function returns 2 components. The first one as explained in the
description is a vector of same length as the first input vector in the
bStrings argument. The second component is a vector that keeps track 
of the time point at which each bit was optimised (see example).
}
\author{
T. Cokelaer
}
\examples{

# Considering the optimised bitstrings at T1, T2 and T3 to be c(1,1,0,1,0,0),
# c(0,1,0) and c(0,1), then we can build the overall bitStrings as follows:

res = buildBitString(list(c(1,1,0,1,0,0), c(0,1,0), c(0,1)))

# The results bit string is accessed through the bs field:
res$bs
#[1] 1 1 0 1 1 1

# and times at which each bit is activated with the bsTimes field:
res$bsTimes

# res$bsTimes = c(1,1,0,1,2,3)

}

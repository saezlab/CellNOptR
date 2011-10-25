\name{normaliseCNOlist}
\alias{normaliseCNOlist}
\title{
Normalisation for boolean modelling.
}
\description{
This function takes in a CNOlist and does the normalisation of the data between 0 and 1, according to two different procedures (see details)
}
\usage{
normaliseCNOlist(CNOlist, EC50Data = 0.5, HillCoef = 2, EC50Noise = 0.1, Detection = 0, Saturation = Inf, ChangeTh = 0, Norm2TorCtrl = "time")
}
\arguments{
  \item{CNOlist}{
a CNOlist
}
  \item{EC50Data}{
parameter for the scaling of the data between 0 and 1, default=0.5
}
  \item{HillCoef}{
Hill coefficient for the scaling of the data, defat to 2
}
  \item{EC50Noise}{
parameter for the computation of a penalty for data comparatively smaller than other time points or conditions
}
  \item{Detection}{
minimum detection level of the instrument, everything smaller will be treated as noise (NA), default to 0
}
  \item{Saturation}{
saturation level of the instrument, everything over this will be treated as NA, default to Inf.
}
  \item{ChangeTh}{
threshold for relative change considered significant, default to 0
}
  \item{Norm2TorCtrl}{
"time" or "ctrl": choice of a normalisation method: compute the relative change compared to the control at the same time, or to the same condition and measurement at time 0
}
}
\details{
The normalisation procedure works as follows: a) every value that is out of the dynamic range of the equipment (as specified by the parameters \code{Detection} and \code{Saturation} are set to NA, b) values are transformed to fold changes relative to the same condition at t0 (if \code{Norm2TorCtrl="time"}) or the control condition (i.e. same inhibitors, no stimuli) at the same time (if \code{Norm2TorCtrl="ctrl"}), c) the fold changes are transformed with a Hill function (i.e. for each data point \code{x^HillCoef/((EC50Data^HillCoef)+(x^HillCoef))}), d) a penalty for "noisiness" is computed for each measurement as the value divided by the maximum value for that readout across all conditions and times (excluding values out of the dynamic range), e) the noise penalty is transformed by a saturation function (for each measurement \code{x/(EC50Noise+x)} where \code{x=x/max(x)}), f) the noise penalty and Hilled fold changes are multiplied, g) if the fold change is negative and bigger than \code{ChangeTh}, the resulting product is multiplied by -1, if the fold change is smaller than \code{ChangeTh} (either positive or negative), it is set to 0.  The normalisation procedure applied here is explained in details in Saez-Rodriguez et al. (2009).

As the normalisation procedure works by computing a fold change relative to the same condition at time 0 or the control condition, if the aforementioned conditions have a value of zero (which is not expected with any common biochemical technique), then the fold change calculation will return a "NaN" value.  If this is a problem for your particular case then we would suggest putting a dummy, very low value, instead of the zero, or setting that measurement to "NA" in the MIDAS file.
}
\value{
a normalised CNOlist
}
\references{
J. Saez-Rodriguez, L. G. Alexopoulos, J. Epperlein, R. Samaga, D. A. Lauffenburger, S. Klamt and P. K. Sorger. Discrete logic modeling as a means to link protein signaling networks with functional analysis of
mammalian signal transduction, Molecular Systems Biology, 5:331, 2009.
}
\author{
C. Terfve
}
\seealso{
makeCNOlist
}
\examples{
#Load a CNOlist

data(CNOlistToy,package="CellNOptR")

#Replace the values in the list by random values 
#(for demonstration purposes, when actually using this function you would simply load a non-normalised CNOlist)

CNOlistToy$valueSignals$t0<-matrix(
	data=runif(n=(dim(CNOlistToy$valueSignals$t0)[1]*dim(CNOlistToy$valueSignals$t0)[2]),min=0,max=400),
	nrow=dim(CNOlistToy$valueSignals$t0)[1],
	ncol=dim(CNOlistToy$valueSignals$t0)[2])
CNOlistToy$valueSignals[[2]]<-CNOlistToy$valueSignals[[1]]+matrix(
	data=runif(n=(dim(CNOlistToy$valueSignals$t0)[1]*dim(CNOlistToy$valueSignals$t0)[2]),min=0,max=100),
	nrow=dim(CNOlistToy$valueSignals$t0)[1],
	ncol=dim(CNOlistToy$valueSignals$t0)[2])

CNOlistToyN<-normaliseCNOlist(
	CNOlistToy, 
	EC50Data = 0.5, 
	HillCoef = 2, 
	EC50Noise = 0.1, 
	Detection = 0, 
	Saturation = Inf, 
	ChangeTh = 0, 
	Norm2TorCtrl = "time")
}
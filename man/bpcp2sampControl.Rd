\name{bpcp2sampControl}
\alias{bpcp2sampControl}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Control function for \code{\link{bpcp2samp}}
}
\description{
Call function to change any one of options, and outputs a list with all defaults except argument that you changed. 
}
\usage{
bpcp2sampControl(Delta = 0, stype = "km", eps = 10^-8,
    nmc=10^6, method="mm.mc", seed=391291)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Delta}{
width of grouped confidence intervals, defaults to 0
}
  \item{stype}{
type of survival estimate, either "km" for Kaplan-Meier or "mue" for median unbiased estimator
}
  \item{eps}{
small value to make integration tractable
}
\item{nmc}{number of Monte Carlo replications}
\item{method}{either 'mm.mc' (method of moments for one sample, meld with Monte Carlo) or 'mc.mc' (Monte Carlo for one sample and melding), see details. }
\item{seed}{random number seed, if NULL do not set random number seed}
}
\details{
We set the seed by default, so that the same data set will always give the same results. If you are doing simulations, this setting of the seed will give problems. So use seed=NULL. 

For method='mm.mc' this is shorthand for calculate the method of moments for one sample, and meld together the two sampling using Monte Carlo methodss; however, technically, that is only done when midp=TRUE, if midp=FALSE then the melding uses numeric integration. 
}
\value{
A list containing the 6 arguments.
}
%\references{}
%\author{}
%\note{}
\seealso{
 \code{\link{bpcp2samp}}
}
\examples{
bpcp2sampControl(Delta=1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{misc}
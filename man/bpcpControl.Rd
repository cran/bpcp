\name{bpcpControl}
\alias{bpcpControl}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Inputs for adjusting numerical calculations in \code{bpcp}
}
\description{
Function that returns a list of arguments.
}
\usage{
bpcpControl(midpMMTol = .Machine$double.eps^0.25, 
    seed=49911,
    tolerance=.Machine$double.eps^0.5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{midpMMTol}{
value used for tol argument in \code{uniroot} call for calculating the midp method of moments method.
}
 \item{seed}{seed for set.seed() when using Monte Carlo method.
If is.null(seed) then do not set the seed.
}
\item{tolerance}{lowest positive value, such that abs(x-y)<tolerance treats x as equal to y. Used in bpcp for seeing if  difference between times are equal to Delta or not.}
}
\details{ 
When doing simulations on the Monte Carlo method, set seed=NULL. Then the seed will not be set at each replication. The default is to set the seed to 49911, so two analyses of the same data on the same version of R will give identical results.

}
\value{
  A list with components named as the arguments.
}
%\references{}
%\author%%  ~~who you are~~}
%\note{
\keyword{optimize}

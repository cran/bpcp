\name{fixtdiff}
\alias{fixtdiff}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Two sample test for Difference in Survival at Fixed Time
}
\description{
Asymptotic two sample tests for difference in survival at a fixed time, using normal approximations and transformations.  See Klien, et al (2007) for details.
}
\usage{
fixtdiff(time,status,group, testtime, 
    trans=c("identity","cloglog","log"),
    varpooled=TRUE, correct=FALSE, doall=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{time}{
time to event for each observation
}
  \item{status}{
status of event time, 1 is observed, 0 is right censored
}
  \item{group}{
group for test, should have two levels, to change order use as factor and change order of levels
}
  \item{testtime}{
fixed time when you want to test for a difference
}
  \item{trans}{type of transformation, one of 'identity', 'cloglog' or 'log'}
  \item{varpooled}{logical, pool the variance?}
\item{correct}{logical, do continuity correction? Continuity correction for when trans='identity' and varpooled (see Warning) }
\item{doall}{logical, do all transformations and corrections}
}
\details{
This function provides p-values for the two sample tests that the survival distributions are equal at time \code{testtime}. The tests are asymptotically normal tests and are described in Klein, et al (2007). These functions are mostly for simulations to evaluate the melded BPCP tests, see \code{\link{bpcp2samp}} and Fay et al (2015).  
}
\value{
A list with the following components:
\item{plo}{one-sided p-value, alternative: S1(testtime)>S2(testtime)}
\item{phi}{one-sided p-value, alternative: S1(testtime)<S2(testtime}
\item{p2}{two-sided p-value, min(1,2*plo,2*phi)}
}
\references{
Fay, MP, Proschan, MA, and Brittain, E (2015) Combining One Sample Confidence Procedures for Inferences in the Two Sample Case. Biometrics 71:146-156. 

Fleiss, Levin, Paik (2003) Statistical Methods for Rates and Proportions, 3rd edition.

Klein, Logan, Harhoff, and Andersen (2007). Analyzing survival curves at a fixed point in time. Statistics in Medicine 26(24): 4505-4519.
}
\author{
Michael P. Fay
}
%\note{}
\section{Warning}{
 Continuity correction derived from the case with no censoriing (see Fleiss et al 3rd edition, pp. 50-55).  May not make sense when there is censoring. Use at own risk. 
}
%\seealso{%% ~~objects to See Also as \code{\link{help}}, ~~~}
\examples{
data(leuk2)
# Note that since the Kaplan-Meier survival at time=35 goes to 
# zero for one group, the results for the log and cloglog 
# transformations are undefined
fixtdiff(leuk2$time,leuk2$status,leuk2$treatment,35,doall=TRUE)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

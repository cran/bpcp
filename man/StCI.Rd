\name{StCI}
\alias{StCI}
\alias{StCI.default}
\alias{StCI.kmci}
\alias{StCI.kmciLR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Get survival and confidence interval at t from \code{kmci}, \code{kmciLR}, or \code{survfit} object}
\description{Just picks out the survival function and confidence interval in a different way depending on the type of object.}
\usage{

\method{StCI}{default}(x,tstar, afterMax = "continue", \dots)


\method{StCI}{kmciLR}(x,tstar, \dots)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a kmci or kmciLR object}
  \item{tstar}{a vector of times that you want survival and CI values}
  \item{afterMax}{character, what to do after tmax (see details)}
\item{\dots}{further arguments to be passed to or from methods.}
}
\details{Since the Kaplan-Meier estimator is undefined after the last observation if it is censored and many confidence interval methods are not defined there either, 
we need to explicitly define what to do. The \code{afterMax} has three options for this: 'continue' (keep surv and ci values the same as the last calculated one), 
'zero' (surv and lower go to zero, upper stays same), 
'half' (surv goes to half value, lower goes to zero, upper stays same).
}
\value{
	The function \code{StCI} returns a data frame with the following variables.
\item{time}{this is tstar}
\item{survival}{survival at tstar}
\item{lower 95 CL}{lower limit at tstar, name changes if different level}
\item{upper 95 CL}{upper limit at tstar, name changes if different level}
}
%\references{  }
\author{Michael Fay
}
%\note{}
%\section{Warning }{} 
\seealso{  \code{\link{kmci}}, \code{\link{kmciLR}} }
\examples{
data(leuk)
## compare to table 1 of Thomas and Grunkmeier (1975)
StCI(kmciTG(leuk$time,leuk$status),c(10,20))
}
\keyword{htest }% at least one, from doc/KEYWORDS
\keyword{survival }% __ONLY ONE__ keyword per line
\keyword{nonparametric}

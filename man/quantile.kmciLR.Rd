\name{quantile.kmciLR}
\alias{quantile.kmciLR}
\alias{quantile.kmci}
\alias{quantile.kmciLRtidy}
\alias{quantile.kmciLRgroup}
\alias{median.kmciLR}
\alias{median.kmci}
\alias{median.kmciLRtidy}
\alias{median.kmciLRgroup}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Quantiles or Medians from \code{kmci}, \code{kmciLR}, \code{kmciLRtidy}, or \code{kmciLRgroup} objects.
}
\description{
Get quantiles or median with the associated confidence intervals from a \code{kmci}, \code{kmciLR}, \code{kmciLRtidy}, or \code{kmciLRgroup} object.
}
\usage{
\method{quantile}{kmciLR}(x, probs = c(0.25, 0.5, 0.75), \dots)
\method{quantile}{kmci}(x, probs = c(0.25, 0.5, 0.75), \dots)
\method{quantile}{kmciLRtidy}(x, probs = c(0.25, 0.5, 0.75), \dots)
\method{quantile}{kmciLRgroup}(x, probs = c(0.25, 0.5, 0.75), \dots)
\method{median}{kmciLR}(x, \dots)
\method{median}{kmci}(x, \dots)
\method{median}{kmciLRtidy}(x, \dots)
\method{median}{kmciLRgroup}(x, \dots)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
a \code{kmci}, \code{kmciLR}, \code{kmciLRtidy}, or \code{kmciLRgroup} object
}
  \item{probs}{
vector of probability to calculate quantiles 
}
  \item{\dots}{
parameters passed
}
}
%\details{}
\value{
A \code{kmciLRtidy} or \code{kmciLRgroup} object will produce a list of matricies. 

The matrix has same number of rows as probs and 4 columns
\item{S(q)}{probs, survival estimate at quantile}
\item{q}{quantile}
\item{lower}{lower confidence limit of q}
\item{upper}{upper confidence limit of q}
}
%\references{%% ~put references to the literature/web site here ~}
%\author{%%  ~~who you are~~}
%\note{%%  ~~further notes~~}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{%% ~~objects to See Also as \code{\link{help}}, ~~~}
\examples{
data(leuk)
data(leuk2)
## kmciTG creates kmci object
fitTG<-kmciTG(leuk$time,leuk$status)
quantile(fitTG)
## bpcp creates kmciLR object 
fitBP<-bpcp(leuk$time,leuk$status)
median(fitBP)

## kmciLRtidy and kmciLRgroups from bpcpfit
practice <- bpcpfit(Surv(time, status)~treatment, data=leuk2)
quantile(practice)
median(practice)

quantile(bpcpfit(leuk2$time, leuk2$status, leuk2$treatment, plotstyle = "standard"))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{survival}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

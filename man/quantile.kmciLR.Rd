\name{quantile.kmciLR}
\alias{quantile.kmciLR}
\alias{quantile.kmci}
\alias{median.kmciLR}
\alias{median.kmci}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Quantiles or Medians from \code{kmci} or \code{kmciLR} objects.
}
\description{
Get quantiles or median with the associated confidence intervals from a \code{kmci} or \code{kmciLR} object.
}
\usage{
\method{quantile}{kmciLR}(x, probs = c(0.25, 0.5, 0.75), \dots)
\method{quantile}{kmci}(x, probs = c(0.25, 0.5, 0.75), \dots)
\method{median}{kmciLR}(x, \dots)
\method{median}{kmci}(x, \dots)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
a \code{kmci} or \code{kmciLR} object
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
matrix same number of rows as probs and with  4 columns
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
## kmciTG creates kmci object
fitTG<-kmciTG(leuk$time,leuk$status)
quantile(fitTG)
## bpcp creates kmciLR object 
fitBP<-bpcp(leuk$time,leuk$status)
median(fitBP)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{survival}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

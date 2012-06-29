\name{summary.kmciLR}
\alias{summary.kmciLR}
\alias{summary.kmci}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summay method for \code{kmci} or \code{kmciLR} object.
}
\description{
Creates a data frame with time (for kmci) or time interval (for kmciLR), survival, lower and upper pointwise confidence intervals. 
}
\usage{
\method{summary}{kmciLR}(object, ...)
\method{summary}{kmci}(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
kmci or kmciLR object
}
  \item{\dots}{
extra arguments
}
}
%\details{}
\value{
creates a data frame. See description.
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
summary(fitTG)
## bpcp creates kmciLR object 
fitBP<-bpcp(leuk$time,leuk$status)
summary(fitBP)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{survival }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

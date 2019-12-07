\name{summary.kmciLR}
\alias{summary.kmciLR}
\alias{summary.kmci}
\alias{summary.kmciLRtidy}
\alias{summary.kmciLRgroup}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summary method for \code{kmci}, \code{kmciLR}, \code{kmciLRtidy}, or \code{kmciLRgroup} object.
}
\description{
Creates a data frame with time (for kmci) or time interval (for kmciLR, kmciLRtidy, and kmciLRgroup), survival, lower and upper pointwise confidence intervals. For kmciLRtidy and kmciLRgroup objects, the group (treatment) variable is also included. 
}
\usage{
\method{summary}{kmciLR}(object, ...)
\method{summary}{kmci}(object, ...)
\method{summary}{kmciLRtidy}(object, ...)
\method{summary}{kmciLRgroup}(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
kmci, kmciLR, kmciLRtidy, or kmciLRgroup object
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
data(leuk2)
## bpcpfit creates kmciLRtidy or kmciLRgroup object
fitBPsurv <- bpcpfit(Surv(time, status)~treatment, data=leuk2)
summary(fitBPsurv)
summary(Surv(time, status)~treatment, data=leuk2, plotstyle = "standard")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{survival}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

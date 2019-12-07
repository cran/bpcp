\name{print.kmciLRtidy}
\alias{print.kmciLRtidy}
\alias{print.kmciLRgroup}
\alias{print.kmciLR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
  Print A Short Summary of a kmciLRtidy, kmciLRgroup, or kmciLR Object
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Print number of observations, number of events, and the median survival with confidence limits for the median of a kmciLRtidy, kmciLRgroup, or kmciLR object with pointwise confidence intervals. The confidence limits will match those specified in the bpcp or bpcpfit function.
}
\usage{
\method{print}{kmciLR}(x, ...)
\method{print}{kmciLRtidy}(x, ...)
\method{print}{kmciLRgroup}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
a kmciLRtidy, kmciLRgroup, or kmciLR obejct
}
  \item{\dots}{
extra arguments
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~}
%\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...}
%\references{
%% ~put references to the literature/web site here ~}
%\author{
%%  ~~who you are~~}
%\note{
%%  ~~further notes~~}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~}
\examples{
library(bpcp)
data(leuk2)
practice <- bpcpfit(Surv(time, status)~treatment, data=leuk2)
practice

bpcpfit(Surv(time, status)~1, data=leuk2)

bpcpfit(Surv(time, status)~1, data=leuk2, plotstyle = "standard")
bpcpfit(leuk2$time, leuk2$status, leuk2$treatment, plotstyle = "standard")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{survival}% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

\name{bpcp-package}
\alias{bpcp-package}
\docType{package}
\title{
Pointwise Confidence Intervals Associated with the Survival Distribution for right censored data
}
\description{
These functions give several different methods for calculating pointwise confidence intervals for the survival distribution for right censored data. 
The primary purpose of the package is to show that the beta product confidence procedure  (using \code{\link{bpcp}}) has better coverage than the best alternatives. 
}
\details{
\tabular{ll}{
Package: \tab bpcp\cr
Type: \tab Package\cr
Version: \tab 0.5.0\cr
Date: \tab 2012-06-28\cr
License: \tab GPL2\cr
LazyLoad: \tab yes\cr
}
The most important function is the  \code{\link{bpcp}} which gives confidence intervals for the survival distribution for right censored data with nice small sample 
properties. That function creates an \code{kmciLR} object which has 6 methods: \code{summary} (create a data frame with estimates and confidence intervals), 
\code{plot} (plot Kaplan-Meier with confidence intervals), \code{lines} (add confidence intervals to a plot),  
\code{StCI} (pick out survival and confidence interval at specific times),
\code{median} (pick out median and confidence interval on median), and \code{quantile} (pick out any quantile and its confidence interval). 
}
\author{
Michael Fay

Maintainer: Michael Fay <mfay@niaid.nih.gov>

}
\references{
Fay, MP, Brittain, E, and Proschan, MA. Pointwise Confidence Intervals for a Survival Distribution for Right Censored Data with Small Samples or Heavy Censoring (unpublished manuscript) 
}
\keyword{ package }
\keyword{survival}
\keyword{nonparametric}
\seealso{
 \code{\link{bpcp}} }
\examples{
data(leuk)
## since there are ties at time=6
## and the data are truncated to the nearest integer, use Delta=1
bfit<-bpcp(leuk$time,leuk$status,Delta=1)
## plot Kaplan-Meier and 95 pct Beta Product Confidence Intervals 
plot(bfit,xlab="time (weeks)")
## details
summary(bfit)
quantile(bfit)
StCI(bfit,2)
}

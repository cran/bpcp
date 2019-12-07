\name{bpcp-package}
\alias{bpcp-package}
\docType{package}
\title{
Pointwise Confidence Intervals Associated with the Survival Distribution for Right Censored Data
}
\description{
The package has functions to give several different methods for calculating pointwise confidence intervals for a single survival distribution for right censored data. There is also a two-sample test for dissimilarity (measured by difference, ratio, or odds ratio) between two survival distributions at a fixed point in time. 


The recommended confidence interval for a single sample is the beta product confidence procedure  (using \code{\link{bpcp}}), and the recommended test for the two-sample test is the melded BPCP test (using \code{\link{bpcp2samp}}).
 
Other confidence intervals and two sample tests are included in the package primarily to compare them to the recommended ones. (And justify the recommendations). 

Also included is a difference in medians test that applies only to non-censored data and is designed to guarantee coverage for all sample sizes (see \code{\link{mdiffmedian.test}}). The test makes no assumptions about the distributions, so that, unlike the Hodges-Lehmann method, tied data are allowed and a shift assumption is not needed.
 }
\details{
\tabular{ll}{
Package: \tab bpcp\cr
Type: \tab Package\cr
Version: \tab 1.3.5\cr
Date: \tab 2017-03-31\cr
License: \tab GPL2\cr
LazyLoad: \tab yes\cr
}
The most important function for the single sample case is the  \code{\link{bpcp}} which gives confidence intervals for the survival distribution for right censored data with nice small sample 
properties. That function creates an \code{kmciLR} object which has 6 methods: \code{summary} (create a data frame with estimates and confidence intervals), 
\code{plot} (plot Kaplan-Meier with confidence intervals), \code{lines} (add confidence intervals to a plot),  
\code{StCI} (pick out survival and confidence interval at specific times),
\code{median} (pick out median and confidence interval on median), and \code{quantile} (pick out any quantile and its confidence interval). A mid-p option for \code{\link{bpcp}} is now available. It gives closer to nominal coverage than the standard (midp=FALSE) BPCP. For details see Fay et al (2013) on the standard BPCP and Fay and Brittain (2016) on the mid-p option. 

For the two-sample test see \code{\link{bpcp2samp}}. This test for equality reduces to Fisher's exact test when there is no censoring. When there is censoring, the test is expected to maintain at least nominal coverage.  For details see Fay et al (2015). 
}
\author{
Michael P. Fay

Maintainer: Michael P. Fay <mfay@niaid.nih.gov>

}
\references{
Fay, MP, Brittain, E, and Proschan, MA. (2013). Pointwise Confidence Intervals for a Survival Distribution with Small Samples or Heavy Censoring.  Biostatistics 14(4): 723-736
doi: 10.1093/biostatistics/kxt016. (copy available at http://www.niaid.nih.gov/about/organization/dcr/brb/staff/Pages/michael.aspx).  

Fay, MP, Proschan, MA, and Brittain, E (2015). Combining One Sample Confidence Procedures for Inference in the Two Sample Case. Biometrics 71:146-156.

Fay, MP, and Brittain, E (2016). Finite Sample Pointwise Confidence Intervals for a Survival Distribution with Right-Censored Data. Statistics in Medicine. doi: 10.1002/sim.6905.
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

\name{bpcp2samp}
\alias{bpcp2samp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Melded BPCP test
}
\description{
Tests for dissimilarity between two groups in their survival distributions at a fixed point in time. Can operationalize that dissimilarity as 'difference', 'ratio' or 
'oddsratio'.
}
\usage{
bpcp2samp(time, status, group, testtime, 
    parmtype = c("difference", "oddsratio", "ratio"), 
    nullparm = NULL, 
    alternative = c("two.sided", "less", "greater"), 
    conf.level = 0.95, 
    control = bpcp2sampControl(Delta = 0, 
       stype = "km", eps = 10^-8))
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
  \item{parmtype}{
parameter type for comparing the survival function of the two groups, either 'difference' 'ratio' or 'oddsratio'
}
  \item{nullparm}{
null value of the parameter of interest, default of NULL gives 0 if parmtype='difference' and 1 otherwise
}
  \item{alternative}{
character, either 'two.sided','less', or 'greater' 
}
  \item{conf.level}{
confidence level, e.g., 0.95
}
  \item{control}{
list of control parameters, see \code{\link{bpcp2sampControl}}
}
}
\details{
The melded confidence interval method is a very general procedure to create confidence intervals for the two sample tests by combining one sample confidence intervals. If S1 and S2 are the survival value at \code{testtime} from sample 1 (first value of \code{group}) and sample 2 (second value of \code{group}) respectively, then 
we can get confidence intervals on the S2-S1 (\code{parmtype='difference'}), S2/S1 (\code{parmtype='ratio'}), or 
(S2*(1-S1))/(S1*(1-S2)) (\code{parmtype='oddsratio'}). 

The resulting melded CIs appear to guarantee coverage as long as the one sample confidence intervals from which the melded CIs are derived have guaranteed coverage themselves. So since we use the BPCP for the one sample intervals and they appear to guarantee coverage (see Fay, Brittain, and Proschan, 2013), we expect the melded BPCP intervals to have at least nominal coverage.  For details see Fay, Proschan and Brittain (2014).   
}
\value{
A list with class "htest" containing the following components:
\item{statistic}{estimate of S1, survival at testtime for group 1}
\item{parameter}{estimate of S2, survival at testtime for group 2}
\item{p.value}{p-value for the test}
\item{conf.int}{a confidence interval for the parameter determined by \code{parmtype}}
\item{estimate}{estimate of parameter determined by \code{parmtype}}
\item{null.value}{the specificed null hypothesized value of the parameter determined by \code{parmtype}}
\item{alternative}{type of alternative with respect to the null.value, either 'two.sided', 'greater' or 'less'}
\item{method}{a character string describing the test}
\item{data.name}{a character string describing the parameter determined by \code{parmtype}}
}
\references{
Fay, MP, Brittain, E, and Proschan, MA. (2013). Pointwise Confidence Intervals for a Survival Distribution with Small Samples or Heavy Censoring.  Biostatistics 14(4): 723-736
doi: 10.1093/biostatistics/kxt016. (copy available at http://www.niaid.nih.gov/about/organization/dcr/brb/staff/Pages/michael.aspx).  



Fay, MP, Proschan, MA, and Brittain, E (2014) Combining One Sample Confidence Procedures for Inferences in the Two Sample Case (unpublished manuscript). 
}
\author{
Michael P. Fay
}
\examples{
data(leuk2)
bpcp2samp(leuk2$time,leuk2$status,leuk2$treatment,35,parmtype="ratio")

bpcp2samp(leuk2$time,leuk2$status,leuk2$treatment,35,parmtype="difference")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

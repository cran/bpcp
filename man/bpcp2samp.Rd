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
    midp=FALSE,
    control = bpcp2sampControl())
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
\item{midp}{logical, do mid-p tests and confidence intervals?}
  \item{control}{
list of control parameters, see \code{\link{bpcp2sampControl}}
}
}
\details{
The melded confidence interval method is a very general procedure to create confidence intervals for the two sample tests by combining one sample confidence intervals. If S1 and S2 are the survival value at \code{testtime} from sample 1 (first value of \code{group}) and sample 2 (second value of \code{group}) respectively, then 
we can get confidence intervals on the S2-S1 (\code{parmtype='difference'}), S2/S1 (\code{parmtype='ratio'}), or 
(S2*(1-S1))/(S1*(1-S2)) (\code{parmtype='oddsratio'}). 

The resulting melded CIs appear to guarantee coverage as long as the one sample confidence intervals from which the melded CIs are derived have guaranteed coverage themselves. So since we use the BPCP for the one sample intervals and they appear to guarantee coverage (see Fay, Brittain, and Proschan, 2013), we expect the melded BPCP intervals to have at least nominal coverage. Note that when there is no censoring the melded CIs derived from the one-sample BPCPs, give matching inferences to Fisher's exact test (i.e., give theoretically identical p-values) when testing the null hypothesis of equality (S1=S2).  For details see Fay, Proschan and Brittain (2015).   

The original melded CIs focused on combining one sample CIs that that guarantee coverage. We can apply the melding to other CIs as well, such as the mid-p style CIs. The mid-p CIs are not designed to guarantee coverage, but are designed to have close to the nominal coverage 'on average' over all the possible values of the parameters. The usual p-value is derived from Pr[ see observed data or more extreme under null], while the mid p-value version comes from (1/2) Pr[see obs data] + Pr[ see more extreme data]. Mid-p CIs come from inverting the test that uses the mid p-value instead of the usual p-value. 



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



Fay, MP, Proschan, MA, and Brittain, E (2015) Combining One Sample Confidence Procedures for Inferences in the Two Sample Case. Biometrics 71:146-156. 
}
\author{
Michael P. Fay
}
\examples{
data(leuk2)
# test ratio of S(20) in two treatment groups, 
# S(20)=survival function at 20 weeks 
bpcp2samp(leuk2$time,leuk2$status,leuk2$treatment,20,parmtype="ratio")
# test difference of S(20) values
bpcp2samp(leuk2$time,leuk2$status,leuk2$treatment,20,parmtype="difference")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}
\keyword{survival}% __ONLY ONE__ keyword per line

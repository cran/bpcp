\name{bpcp-internal}
\alias{abmm}
\alias{kmgw.calc}
\alias{kmci.mid}
\alias{kmci.cons}
\alias{kmConstrain}
\alias{kmcilog}
\alias{kmConstrainBeta.calc}
\alias{bpcp.mm}
\alias{bpcp.mc}
\alias{bpcpMidp.mm}
\alias{getmarks}
\alias{getmarks.x}
\alias{intChar}
\alias{citoLR}
\alias{meldMC}
\alias{betaMeldTestMidp.mc}

\alias{borkowf.calc}

\alias{qqbeta}
\alias{uvab}
\alias{rejectFromInt}


%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Internal functions
}
\description{
Functions  called by other functions. Not to be directly called by user.
}
\usage{
abmm(a1,b1,a2,b2)
kmgw.calc(time, status, keepCens = TRUE)
borkowf.calc(x, type = "log", alpha = 0.05)
kmConstrain(tstar, pstar, x, alpha = 0.05)
kmConstrainBeta.calc(tstar, pstar, x, alpha = 0.05)
bpcp.mm(x,alpha=0.05)
bpcp.mc(x,nmc=100,alpha=0.05, testtime=0, DELTA=0, midp=FALSE)
bpcpMidp.mm(x,alpha=0.05, midptol=.Machine$double.eps^0.25)
kmcilog(x, alpha = 0.05)

qqbeta(x, a, b)
rejectFromInt(theta,interval,thetaParm=FALSE)
uvab(u, v)
citoLR(x)

getmarks(time, status)
getmarks.x(x)


intChar(L, R, Lin = rep(FALSE, length(L)), Rin = rep(TRUE, length(L)), digits = NULL)


meldMC(T1,T2, nullparm=NULL, 
    parmtype=c("difference","oddsratio","ratio","cdfratio"),
    conf.level=0.95, 
    alternative=c("two.sided","less","greater"),
    dname="",estimate1=NA, estimate2=NA)

betaMeldTestMidp.mc(betaParms1,
    betaParms2,nullparm=NULL, 
    parmtype=c("difference","oddsratio","ratio","cdfratio"),
    conf.level=0.95, conf.int=TRUE,
    alternative=c("two.sided","less","greater"),
    dname="",
    estimate1=NA, estimate2=NA, nmc=10^6)


}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{a}{ beta shape1 parameter}
  \item{b}{ beta shape2 parameter }
  \item{a1}{first beta shape1 parameter, first of two beta distributions}
  \item{a2}{second beta shape1 parameter, second of two beta distributions}
  \item{b1}{first beta shape2 parameter, first of two beta distributions}
  \item{b2}{second beta shape2 parameter, second of two beta distributions}
  \item{u}{vector of means of beta distributions}
  \item{v}{vector of variances of beta distributions}
  \item{time}{time to event or censoring}
  \item{status}{vector of event status, 1 for events 0 for censoring}
  \item{keepCens}{logical, keep times with only censored values?}
  \item{x}{output from kmgw.calc}
  \item{theta}{either the parameter under the null (if thetaParm=TRUE) or an estimate of theta (if thetaParm=FALSE)}
  \item{thetaParm}{logical, is theta a parameter?}
  \item{interval}{either a confidence interval (if thetaParm=TRUE) or quantiles from a null distribution (if thetaParm=FALSE)}
  \item{alpha}{1-conf.level}
 \item{testtime}{time for test, needed for output for two-sample test}
\item{midp}{logical, do mid-p tests and/or confidence intervals?}
\item{midptol}{tol value passed to uniroot in function}
\item{DELTA}{same at Delta in bpcp}
  \item{tstar}{time for survival distribution}
  \item{pstar}{null value for survival}
  \item{type}{character describing method, either 'log' transformation, 'logs' log transformation with shift, 'norm' no transformation, 'norms' no transformation with shift }
\item{nmc}{number of Monte Carlo reps}
\item{L}{left end of  intervals associated with each surv and ci value}
\item{R}{right end of  intervals associated with each surv and ci value}
\item{Lin}{logical vector, include left end in interval?}
\item{Rin}{logical vector, include right end in  interval?}
\item{digits}{how many significant digits to use}
\item{T1}{vector of nmc simulated values for parameter from group 1}
\item{T2}{vector of nmc simulated values for parameter from group 2}
\item{nullparm}{null value of the 2 sample parameter, when NULL gives values appropriate for parmtype}
\item{parmtype}{type of parameter for the two sample test, for details see \code{\link{bpcp2samp}}  }
\item{conf.level}{confidence level}
\item{conf.int}{logical, calculate confidence interval?}
\item{alternative}{alternative hypothesis}
\item{dname}{data name for 'htest' class of the result}
\item{estimate1}{estimate of parameter from group 1}
\item{estimate2}{estimate of parameter from group 2}
\item{betaParms1}{named list of beta parameters from group 1 (usually come from method of moments), names: alower,blower, aupper, bupper}
\item{betaParms2}{named list of beta parameters from group 2, names: alower,blower, aupper, bupper}
}
\details{
abmm uses method of moments to find a,b parameters from beta distribution that is product of two other beta RVs.   

kmgw.calc calculates the Kaplan-Meier and Greenwood variances.

kmci.mid and kmci.cons calculate confidence intervals using a new method with either mid-p-like intervals or a conservative interval
from input from kmgw.calc.
 
borkowf.calc calculates the Borkowf intervals from output from kmgw.calc.

kmcilog gives normal approximation confidence interval using log transformation.

bpcp.mm and bpcp.mc are the main calculation functions (.mm for method of moments, .mc for Monte Carlo simulation) 
for bpcp (repeated Beta method). Both output a list with two vectors, upper and lower. bpcpMidp.mm and bpcpMidp.mc are the mid-p versions of these functions. 

kmConstrain gives constrained K-M estimate, and kmConstrainBeta.calc gives ci and tests using Beta distribution. 
 
qqbeta is like qbeta, but allows a=0 (giving a value of 0 when b>0) and b=0 (giving a value of 1 when a>0).

rejectFromInt inputs theta and an interval and gives a vector with 3 terms, estGTnull=1 if reject and estimate is greater than null value,
estLTnull=1 if reject and estimate is less than null value, two.sided=1 if reject in either direction. 
The thetaParm=TRUE means that theta is the parameter under the null so that interval is a confidence interval,
while thetaParm=FALSE means that theta is an estimate of the parameter and interval are quantiles from the null distribution. 

uvab takes means and variances of beta distributions and returns shape parameters. 
}
%\value{}
%\references{}
\author{Michael Fay}
%\note{%%  ~~further notes~~}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{%% ~~objects to See Also as \code{\link{help}}, ~~~}
%\examples{}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{internal}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

\name{kmtestALL}
\alias{kmtestALL}
\alias{kmciBorkowf}
\alias{kmtestConstrainBoot}
\alias{kmtestConstrainBeta}
\alias{bpcp}
\alias{kmciSW}
\alias{kmciTG}
\alias{kmci1TG}
\alias{kmtestBoot}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Pointwise confidence intervals for survival for right censored data.
}
\description{
These functions give several different methods for calculating pointwise confidence intervals for the survival distribution for right censored data. 
The recommended confidence intervals are the beta product ones given by \code{bpcp}. 

The other confidence intervals are included primarily to show that the beta product confidence procedure (using \code{bpcp}) has better coverage than the best alternatives. 
See details for a description of all the methods.
}
\usage{
bpcp(time,status,nmc=0,alpha=.05,Delta=0,stype="km")


kmciBorkowf(time, status, type="log", alpha = 0.05)
kmtestBoot(time, status, tstar, pstar, M = 1000, alpha = 0.05)
kmtestConstrainBoot(time, status, tstar, pstar, M = 1000, alpha = 0.05)
kmtestConstrainBeta(time,status,tstar,pstar,alpha=.05)
kmciSW(time, status, alpha = 0.05)
kmciTG(time, status, alpha = 0.05)
kmci1TG(time, status, tstar, alpha = 0.05)

kmtestALL(time,status,t0,S0,cens=NULL,M=1000,NMC=10^5,alpha=0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{time}{
time to event or censoring
}
  \item{status}{ status vector,
1 is event, 0 is censoring
}
  \item{alpha}{
1- conf.level
}
\item{nmc}{number of Monte Carlo replications from each beta distribution, nmc=0 means use method of moments for beta parameters instead}
\item{NMC}{same as nmc}
\item{Delta}{width of grouped confidence intervals, defaults to 0}
\item{stype}{type of survival estimate, either "km" for Kaplan-Meier or "mue" for median unbiased estimator} 
  \item{tstar}{ time to test survival distribution}
  \item{pstar}{null survival distribution}
  \item{M}{number of bootstrap replications}
\item{t0}{null hypothesis time for survival test }
\item{S0}{null hypothesis value of survival at t0}
\item{cens}{vector of censoring times (even those with failures before it), used for Binomial test. If NULL gives NA for binom test}
  \item{type}{see details}

}
\details{ 

The beta product confidence procedure will give pointwise confidence intervals for right censored data with the following properties. When there is no censoring 
or Progressive Type II censoring the BPCP guarantees central coverage (e.g., the error rate on either side of the 95 percent confidence interval is guaranteed to be less than 2.5 percent). 
For general independent censoring the BPCP is asymptotically equivalent to standard methods such as the normal approximation with Greenwood variance, and 
hence the BPCP (as with the other confidence interval given here) goes to the correct confidence interval for any t.   

Now we describe the other methods. 
In general the functions are of three naming types: kmtestXX, kmci1XX and kmciXX, where XX changes for different methods. 
Functions kmtestXX only test whether S(tstar)=pstar and return a vector of 1s for reject and 0s for fail to rejecting either of the one-sided or the two-sided hypotheses. 
Functions kmci1XX only give confidence intervals at  S(tstar), while kmciXX give confidence intervals for all values of t.
The standard methods calculate the confidence intervals at the observed failure times and carry them forward (e.g., kmciTG, kmciSW)
and the results are objects of class kmci. More involved methods allow confidence intervals to change after censored objects (e.g., kmciBorkowf, bpcp)
and the results are objects of class kmciLR.  

The function kmtestBoot tests S(tstar)=pstar using the nonparametric bootstrap (sampling vectors of (time,status) with replacement) with the percentile method as described in Efron (1981). 
The function kmtestConstrainBoot and kmtestConstrainBeta tests S(tstar)=pstar using the constrained Bootstrap or constrained Beta method described 
in Barber and Jennison (1999). 

The function kmci1TG does a confidence interval only at tstar, while kmciTG does a confidence interval 
at all the observed event times. The method can be derived as a likelihood ratio test and is described in  Thomas and Grunkemeier (1975). 
It has asymptoticly correct coverage, which is  rigorously proved in 
Murphy (1995). You can also think of the method as the empirical likelihood applied to the survival distribution for 
right censored  data (see Owen, 2001, p. 144-145).


The function \code{kmciSW} calculates confidence intervals using Edgeworth expansions as described in 
Strawderman and Wells (1997). Note, Strawderman, Parzen and Wells (1997) is easier to understand than  Strawderman and Wells (1997). 


Borkowf (2005) creates confidence intervals for the Kaplan-Meier survival estimate for right censored data. 
He allows the confidence interval to change at censoring times as well as at failure times.  

Four types of confidence intervals may be selected. The asymptotic normal approximation (type="norm"), the shifted K-M 
estimate with normal approximation (type="norms"), the log transformed normal approximation using the delta method (type="log"), 
and the log transformed normal approximation using the delta method with the shifted K-M  (type="logs"). 


The function \code{kmtestALL} performs hypothesis tests on all the methods except the unconstrained bootstrap method (unless M=0 then it does not test the constrained bootstrap method either). 
The output is a matrix with three columns with a value of 1 representing either (1) rejection for two-sided test implying the estimate is greater than the null, 
(2) rejection for two-sided test implying the estimate less than the null, or (3) any rejection of the two-sided test. 
Each row represents a different test.


The \code{kmci} and \code{kmciLR} class have the following methods: \code{summary}, \code{plot}, \code{lines}, \code{quantile}, and \code{median}.
}
\value{
The functions return an object of class either kmci or kmciLR (see details). 

both kmci and kmciLR objects are lists, and both contain elements
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{surv}{survival distribution in interval}
\item{lower}{lower pointwise confidence limit in interval}
\item{upper}{upper pointwise confidence limit in interval}

additionally the kmci objects have an element
\item{time}{time of survival or confidence interval}
while the kmciLR have intervals represented by the four elements
\item{L}{left endpoint of interval}
\item{Lin}{logical vector, include left endpoint?}
\item{R}{right endpoint of interval}
\item{Rin}{logical vector, include right endpoint?}
}
\references{
Fay, MP, Brittain, E, and Proschan, MA. (2013). Pointwise Confidence Intervals for a Survival Distribution with Small Samples or Heavy Censoring.  Biostatistics (to appear)
doi: 10.1093/biostatistics/kxt016. (copy available at http://www.niaid.nih.gov/about/organization/dcr/brb/staff/Pages/michael.aspx).  


Barber and Jennison (1999) Biometrics, 55: 430-436.

Borkowf (2005) Statistics in Medicine, 24: 827-851.

Efron (1981) JASA 76:312-319.

Murphy (1995) JASA 90: 1399-1405.

Owen (2001) Empirical Likelihood. Chapman and Hall: New York.

Strawderman and Wells (1997) JASA 92:1356-1374.

Strawderman, Parzen and Wells (1997) Biometrics 53: 1399-1415.

Thomas and Grunkemeier (1975) JASA 70: 865-871.


}
\author{
Michael Fay
}
%\note{%%  ~~further notes~~}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{%% ~~objects to See Also as \code{\link{help}}, ~~~}
\examples{
library(bpcp)
data(leuk)

### Recommended method is bpcp
### since the data are truncated to the nearest integer
### use Delta=1 option
out<-bpcp(leuk$time,leuk$status,Delta=1)
summary(out)
median(out)
plot(out)



### Borkowf 2005 method
norm<-kmciBorkowf(leuk$time,leuk$status,type="norm")
norms<-kmciBorkowf(leuk$time,leuk$status,type="norms")
## check Table VII of Borkowf
I<-c(1,2,3,5,7,8,9,11,13,15,17,19,21,23,25,27,29,31,33)
round(data.frame(lowerNorm=norm$lower[I],
    upperNorm=norm$upper[I],lowerNormS=norms$lower[I],
    upperNorms=norms$upper[I],row.names=norm$Interval[I]),3)

### Strawderman and Wells (1997) method
swci<-kmciSW(leuk$time,leuk$status)
summary(swci)


### Thomas and Grunkemeier 1975 method
x<-kmciTG(leuk$time,leuk$status)
summary(x)
## compare to Table 1, Sample 2, of Thomas and Grunkemeier (1975)
StCI(x,c(10,20))


}
\keyword{htest }% at least one, from doc/KEYWORDS
\keyword{survival }% __ONLY ONE__ keyword per line
\keyword{nonparametric}

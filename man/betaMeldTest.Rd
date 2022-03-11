\name{betaMeldTest}
\alias{betaMeldTest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Melded Beta Test
}
\description{
This function gives a two sample melded beta test together with the associated melded confidence intervals. It can be used when the confidence distributions (or upper and lower confidence distributions) for the one sample parameters are beta, and one is interested in either the difference, odds ratio, or ratio of those two one sample parameters. The betaMeldTest is usually called by \code{\link{bpcp2samp}}, and not called directly by the user.
}
\usage{
betaMeldTest(betaParms1, betaParms2,
    nullparm = NULL, 
    parmtype = c("difference", "oddsratio", "ratio","cdfratio"), 
    conf.level = 0.95, conf.int = TRUE, 
    alternative = c("two.sided", "less", "greater"), 
    eps = 10^-8, dname = "", estimate1 = NA, estimate2 = NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{betaParms1}{
a list of the beta parameters for group 1
}
  \item{betaParms2}{
a list of the beta parameters for group 2
}
  \item{nullparm}{
null value of the parameter of interest, default of NULL gives 0 if parmtype='difference' and 1 otherwise
}
  \item{parmtype}{
parameter type for comparing the survival function of the two groups, either 'difference' 'ratio' 'oddsratio' or 'cdfratio'
}
  \item{conf.level}{
confidence level, e.g., 0.95
}
  \item{conf.int}{
logical, calculate confidence interval?
}
  \item{alternative}{
character, either 'two.sided','less', or 'greater' 
}
  \item{eps}{
small value to make integration tractable
}
  \item{dname}{
name describing data
}
  \item{estimate1}{
estimate of mean for beta parameter of group 1 (statistic of htest object)
}
  \item{estimate2}{
estimate of mean for beta parameter for group 2 (parameter of htest object)
}
}
\details{
If the upper and lower confidence distributions for both samples are described by beta distributions, then you can create a CD test using this function. For example,
if you have sample 1 is binomial with x (with 0<x<n) out of n positive responses, then the 100(1-alpha) confidence interval is qbeta(alpha/2, x, n-x+1) and 
qbeta(1-alpha/2, x+1, n-x). So the lower confidence distribution is beta with parameters a=x and b=n-x+1, and the upper CD is beta with parameters a=x+1 and b=n-x. 

See \code{\link{bpcp2samp}} for a full description of the parmtypes.
}
\value{
an object of class 'htest'
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\examples{
fisher.test(matrix(c(4,5,2,22),2,2),alternative="greater")

betaMeldTest(
    betaParms1=list(alower=2,blower=22+1,aupper=2+1,bupper=22),
    betaParms2=list(alower=4,blower=5+1,aupper=4+1,bupper=5),
    alternative="greater",parmtype="oddsratio",
    estimate1=2/24,estimate2=4/9)

}
\keyword{htest}

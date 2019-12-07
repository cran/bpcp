\name{mdiffmedian.test}
\alias{mdiffmedian.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Melded Difference in Medians Test
}
\description{
Tests for a difference in two medians. No assumptions about the two distributions are needed (may be discrete with ties allowed, no shift assumption is required). Uses the melded confidence interval derived from the one sample confidence intervals associated with the sign test (a version that allows for ties). Derivation of the test does not require large samples, and confidence intervals are intended to guarantee coverage regardless of sample size. 
}
\usage{
mdiffmedian.test(x1, x2, nulldiff = 0, 
    alternative = c("two.sided", "less", "greater"), 
    conf.level = 0.95)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x1}{
vector of numeric responses from group 1
}
  \item{x2}{
vector of numeric responses from group 2
}
  \item{nulldiff}{
difference in medians under the null, median(x2)-median(x1)
}
  \item{alternative}{
a character string specifying the alternative hypothesis, must be one of "two.sided" (default), "greater" or "less". You can specify just the initial letter.
}
  \item{conf.level}{
confidence level of the interval.
}
}
\details{
The melded confidence interval method is a general method for combining two one-sample confidence intervals (CIs). In this function, we use the melded CI method on the two one-sample CIs from the sign test that allows for ties. This creates CIs for the difference in medians that requires very few assumptions. In particular, ties are allowed and no shift assumption is needed.  For details see Fay, Proschan and Brittain (2015). 
}
\value{
a list of class 'htest' with elements:
\item{statistic}{median of x1}
\item{parameter}{median of x2}
\item{p.value}{p-value of the test}
\item{conf.int}{confidence interval for the difference in medians}
\item{estimate}{median(x2)-median(x1)}
\item{null.value}{null hypothesis value for difference in medians}
\item{alternative}{type of alternative hypothesis}
\item{method}{description of test}
\item{data.name}{description of input}
}
\references{
Fay, MP, Proschan, MA, Brittain, E (2015). Combining One-sample confidence procedures for inference in the two-sample case. Biometrics. 71: 146-156.
}
\author{
Michael P. Fay
}
\note{This function does not allow censoring. Also, there is a price for not needing large samples nor assumptions about the distributions: if you do not have enough data, your confidence intervals may be the entire real line. For example, if you have continuous data with equal sample sizes in both groups, then if you have 6 or fewer observations in each group, then the 95 percent confidence interval on the difference in medians will be (-Inf,Inf).  }
\examples{
set.seed(1)
trtA<-rpois(20,1.5)
trtB<-rpois(23,5.5)
mdiffmedian.test(trtA,trtB)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}
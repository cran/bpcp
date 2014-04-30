\name{leuk2}
\alias{leuk2}
\docType{data}
\title{
Acute Leukemia data from Freireich et al (1963).
}
\description{
 In this study there were 21 pairs of subjects, and within each pair one subject received 6-mercaptopurine (6-MP) and one got placebo. The data are right censored. 

See also Gehan (1965) who used the data ignoring the pairing so that he could illustrate his famous two-sample (non-paired) rank test.
}
\usage{data(leuk2)}
\format{
  A data frame with 42 observations on the following variables.
  \describe{
    \item{\code{time}}{time in remission (in weeks)}
    \item{\code{status}}{event status, 1 is relapse, 0 is censored}
    \item{\code{treatment}}{treatment group: either 'placebo' or '6-MP'}
    \item{\code{pair}}{pair id number}

  }
}
%\details{}
%\source{}
\references{
Gehan (1965) Biometrika 52:203-223.

Freireich et al (1963) Blood 21(6):699-716.
}
\seealso{  \code{\link{leuk}} is only the treated group
}
\examples{
data(leuk2)
}
\keyword{datasets}

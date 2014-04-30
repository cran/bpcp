\name{leuk}
\alias{leuk}
\docType{data}
\title{
Acute Leukemia data (treatment only) from Freireich et al (1963).
}
\description{
 This is only the 21 patients who received 6-mercaptopurine (6-MP). There
were 21 patients who got placebo (see \code{\link{leuk2}} for complete data). 

See also Borkowf (2005)
}
\usage{data(leuk)}
\format{
  A data frame with 21 observations on the following 2 variables.
  \describe{
    \item{\code{time}}{time in remission (in weeks)}
    \item{\code{status}}{event status, 1 is relapse, 0 is censored}
  }
}
%\details{}
%\source{}
\references{
Borkowf (2005) Statistics in Medicine, 24: 827-851.

Freireich et al (1963) Blood 21(6):699-716.
}
\seealso{ \code{\link{leuk2}} for complete data.
}
\examples{
data(leuk)
}
\keyword{datasets}

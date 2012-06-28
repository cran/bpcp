\name{leuk}
\alias{leuk}
\docType{data}
\title{
Acute Leukemia data from Freireich et al (1963).
}
\description{
 This is only the 21 patients who received 6-mercaptopurine (6-MP). There
were 21 patients who got placebo, but their remission time is not given. See also Borkowf (2005)
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
\examples{
data(leuk)
## maybe str(luek) ; plot(luek) ...
}
\keyword{datasets}

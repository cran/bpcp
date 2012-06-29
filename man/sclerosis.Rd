\name{sclerosis}
\alias{sclerosis}
\docType{data}
\title{
Pilot study of treatment of severe systemic sclerosis (Nash, et al, 2007).
}
\description{
Severe systemic sclerosis  is a serious autoimmune disease affecting multiple organs including the heart, lungs, kidney, and skin.  
Between 1997 and 2005, a cohort of 34 patients was enrolled in a single arm pilot study of high-dose immunosuppressive therapy and autologous hetapoietic cell  transplantation
}
\usage{data(sclerosis)}
\format{
  A data frame with 34 observations on the following 3 variables.
  \describe{
    \item{\code{day}}{time to relapse or censoring, in days}
    \item{\code{year}}{time to relapse or censoring, in years (day/365.25)}
    \item{\code{status}}{0 is censored, 1 is event}
  }
}
%\details{}
%\source{%%  ~~ reference to a publication or URL from which the data were obtained ~~}
\references{
Nash, R.A., McSweeney, P.A., Crofford, L.J., Abidi, M., Chen, C.S., Godwin, J.D., Gooley, T.A., Holmberg, L., Henstorf, G., LeMaistre, C.F., others (2007).
"High-dose immunosuppressive therapy and autologous hematopoietic cell transplantation for severe systemic sclerosis: long-term follow-up of the US multicenter pilot study"
Blood 110 (4): 1388-.
}
\examples{
data(sclerosis)
plot(bpcp(sclerosis$year,sclerosis$status))
}
\keyword{datasets}

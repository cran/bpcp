\name{kmci.object}
\alias{kmci.object}
\alias{kmciLR.object}
\alias{kmci}
\alias{kmciLR}

\title{
  Kaplan-Meier (Survival Curve) Confidence Interval Object 
}
\description{
There are two classes of objects described in this help, the \code{kmciLR} class returned by the functions \code{\link{bpcp}} or \code{\link{kmciBorkowf}},
 and the   \code{kmciLR} class returned by the functions  \code{\link{kmciTG}} or \code{\link{kmciSW}}. Both classes 
represent a fitted survival curve with pointwise confidence intervals. 

The \code{kmciLR} class allows for confidence intervals to change at  any time point, while the \code{kmci} class only has the confidence intervals change at observed failures.


Objects of this class has methods for the functions \code{summary}, \code{plot}, \code{lines}. 


Class of objects that represent a fitted survival curve. 
}
\section{Structure}{
The following components must be included in a legitimate 
\code{kmci} 
object. 
}
\arguments{
\item{time}{
the time points of observed failures (assumed surv and lower and upper steps that these times) 
}
\item{cens}{
time points where there is censoring but no observed failure
}
\item{surv}{
the estimate of survival at time t+0. 
This is a vector. 
}
\item{upper}{
upper confidence limit for the survival curve. 
}
\item{lower}{
lower confidence limit for the survival curve. 
}
\item{conf.level}{
the level of the confidence limits, e.g. 90 or 95\%. 
}
}
\seealso{
\code{\link{plot.kmci}},  
\code{\link{summary.kmci}},
\code{\link{StCI.kmci}},
\code{\link{median.kmci}},
\code{\link{quantile.kmci}}.   
}
\keyword{survival}

\name{kmci.object}
\alias{kmci.object}
\alias{kmci}

\title{
  Kaplan-Meier (Survival Curve) Confidence Interval Object 
}
\description{
The \code{kmci} class is returned by the functions  \code{\link{kmciTG}} or \code{\link{kmciSW}}. The  class 
represents a fitted survival curve with pointwise confidence intervals. 

Unlike the \code{\link{kmciLR}} class, which  allows for confidence intervals to change at  any time point, 
the \code{kmci} class only has the confidence intervals change at observed failures.


Objects of this class has methods for the functions \code{summary}, \code{plot}, \code{lines}. 


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
the level of the confidence limits, e.g., 0.95. 
}
}
\seealso{
\code{\link{kmciLR.object}}
\code{\link{plot.kmci}},  
\code{\link{summary.kmci}},
\code{\link{StCI.kmci}},
\code{\link{median.kmci}},
\code{\link{quantile.kmci}}.   
}
\keyword{survival}

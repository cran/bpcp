\name{kmciLR.object}
\alias{kmciLR.object}
\alias{kmciLR}

\title{
  Kaplan-Meier (Survival Curve) Confidence Interval LR Tidy Object 
}
\description{

The \code{kmciLR} class allows for confidence intervals to change at  any time point, 
while the \code{\link{kmci}} class only has the confidence intervals change at observed failures.


Objects of this class has methods for the functions \code{summary}, \code{plot}. 

formula, data, nmc=0, alpha=.05, Delta=0, stype="km", midp=FALSE,
                         monotonic=NULL, control=bpcpControl()

}
\section{Structure}{
The following components must be included in a legitimate 
\code{kmciLR} 
object. 
}
\arguments{
\item{cens}{
time points where there is censoring but no observed failure
}
\item{surv}{
the estimate of survival in the interval described by L and R. This is a vector. 
}
\item{upper}{
upper confidence limit for the survival curve in the interval described by L and R. 
}
\item{lower}{
lower confidence limit for the survival curve in the interval described by L and R. 
}
\item{L}{ 
vector of left ends of interval associated with lower and upper
}
\item{Lin}{ 
vector of logicals, should left end of interval be included?
}
\item{R}{ 
vector of right ends of interval associated with lower and upper
}
\item{Rin}{ 
vector of logicals, should right end of interval be included?
}
\item{Interval}{ 
character vector describing intervals
}
\item{stype}{ 
character vector giving type of survival estimate, either 'km' or 'mue'
}
\item{conf.level}{
the level of the confidence limits, e.g., 0.95. 
}
}
\seealso{
\code{\link{plot.kmciLR}},  
\code{\link{summary.kmciLR}},
\code{\link{bpcp}}.
}
\keyword{survival}

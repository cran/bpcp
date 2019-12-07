\name{kmciLRtidy.object}
\alias{kmciLRtidy.object}
\alias{kmciLRtidy}

\title{
  Kaplan-Meier (Survival Curve) Confidence Interval LR Tidy Object 
}
\description{
The \code{kmciLRtidy} class is returned by the function \code{\link{bpcpfit}} when the plotstyle = "ggplot", and
represents a fitted survival curve with pointwise confidence intervals. This object may contain more than one
fitted survival curve with the corresponding pointwise confidence intervals based on a treatment/grouping variable, 
if one is specified in the formula of \code{\link{bpcpfit}}.

The \code{kmciLRtidy} class allows for confidence intervals to change at any time point (similar to the \code{\link{kmciLR}} class), 
while the \code{\link{kmci}} class only has the confidence intervals change at observed failures.

Objects of this class has methods for the functions \code{summary}, \code{plot}, \code{print}, \code{quantile}, and \code{median}.

When printed, objects of this class will display the total number of subjects, total number of events, median survival time, and the confidence limits(for the pre-specified confidence level), by treatment group if applicable.

}
\section{Structure}{
A \code{kmciLRtidy} object is a list of kmciLR objects (also lists). The length of the kmciLRtidy object corresponds to the number of treatment/grouping variables; each item in the list is a kmciLR object for the given treatment (the name of each item in the list is the name of each group). If no grouping variable is given, the output is a list of 1 kmciLR object. 

The following components must be included in each element of a legitimate 
\code{kmciLRtidy} 
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
\item{num}{ 
total number of subjects at the start of that survival curve.
}
\item{events}{ 
total number of events (observed failures) for that survival curve. 
}

}
\seealso{
\code{\link{plot.kmciLRtidy}},  
\code{\link{summary.kmciLRtidy}},
\code{\link{print.kmciLRtidy}}, 
\code{\link{quantile.kmciLRtidy}},
\code{\link{median.kmciLRtidy}},
\code{\link{bpcpfit}}.
}
\keyword{survival}

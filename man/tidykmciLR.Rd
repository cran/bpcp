\name{tidykmciLR}
\alias{tidykmciLR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
  Dataframe of \code{kmciLRtidy}, \code{kmciLRgroup}, or a \code{kmciLR} object.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
  Takes a \code{kmciLRtidy}, \code{kmciLRgroup}, or a \code{kmciLR} object (a list) and converts it into a drameframe, which can further be used in plotting. Every two time points represents a time interval.  
}
\usage{
tidykmciLR(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
a kmciLRtidy, kmciLRgroup, or kmciLR object
}
}
\details{
creates a dataframe. See description. 
}
\value{
This function returns a dataframe with the following columns:

\item{time}{Time}
\item{surv}{Value of survival curve at that time point}
\item{lower}{ Lower bound of the CI for the survival curve}
\item{upper}{ Upper bound of the CI for the survival curve} 
\item{group}{treatment or grouping variable (if applicable)}

There are two rows per time point representing the change in either the survival function or confidence bands.
}

\examples{
library(bpcp)
data(leuk2)
practice <- bpcpfit(Surv(time, status)~treatment, data=leuk2)

tidy <- tidykmciLR(practice)
ggplot(tidy, aes(x = time, y = surv, ymin = lower, ymax = upper, col = group)) + 
  geom_line(show.legend=FALSE) + geom_ribbon(alpha = .2, aes(fill=group)) + xlab("Time") + 
  ylab("Survival") + ggtitle("K-M curves with bpcp CIs") 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{survival}% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

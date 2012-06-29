\name{plot.kmciLR}
\alias{plot.kmci}
\alias{plot.kmciLR}
\alias{lines.kmci}
\alias{lines.kmciLR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot and lines methods for \code{kmci} and \code{kmciLR} objects. 
}
\description{
Plots survival curves and/or confidence intervals. 
}
\usage{

\method{plot}{kmci}(x, ...)

\method{plot}{kmciLR}(x, XLAB = "time", YLAB = "Survival", YLIM = c(0, 1), ciLTY = 2, ciCOL = gray(0.8), mark.time = NULL, linetype = "both", ...)

\method{lines}{kmciLR}(x, lty = c(2, 1), col = c(gray(0.8), gray(1)), linetype = "ci", mark.time = NULL, ...)

\method{lines}{kmci}(x, ...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
\code{kmci} or \code{kmciLR} object (created by functions described in \code{\link{kmtestALL}}) 
}
  \item{XLAB}{
label for x axis
}
  \item{YLAB}{
label for y axis
}
  \item{YLIM}{
limits for y axis
}
  \item{ciLTY}{
lty (line type) for confidence intervals
}
  \item{ciCOL}{
col (color) for confidence intervals
}
\item{col}{
vector of colors, first element used for ci second for survival curve
}
\item{lty}{
vector of line types, first element used for ci second for survival curve
}
  \item{mark.time}{
put hash marks for censored objects (default puts marks of stype="km" but not if stype="mue")
}
  \item{linetype}{
character, which lines to draw: either 'both', 'surv' or 'ci'
}
  \item{\dots}{
Extra parameters to be passed. Any argument in \code{plot.kmciLR} can be passed from \code{plot.kmci}, similarly for \code{line}. 
Other parameters are usually graphical parameters passed to plot and segment calls within function.
}
}
%\details{}
%\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
%}
%\references{
%% ~put references to the literature/web site here ~
%}
%\author{
%%  ~~who you are~~
%}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
data(leuk)
## kmciTG creates kmci object
fitTG<-kmciTG(leuk$time,leuk$status)
plot(fitTG)
## bpcp creates kmciLR object 
fitBP<-bpcp(leuk$time,leuk$status)
lines(fitBP,lwd=3,lty=1,col=gray(.5))
legend(0,.2,legend=c("Kaplan-Meier","Thomas-Grunkemeier 95 pct CI","Beta Product 95 pct CI"),lwd=c(1,1,3),lty=c(1,2,1),col=c(gray(0),gray(.8),gray(.5)))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{survival}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

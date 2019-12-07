\name{plot.kmciLR}
\alias{plot.kmci}
\alias{plot.kmciLR}
\alias{plot.kmciLRtidy}
\alias{plot.kmciLRgroup}
\alias{lines.kmci}
\alias{lines.kmciLR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot and lines methods for \code{kmci}, \code{kmciLR}, \code{kmciLRtidy}, and \code{kmciLRgroup} objects. 
}
\description{
Plots survival curves and/or confidence intervals. \code{kmciLR} objects from the bpcp function will plot using base R plot. 
If a \code{kmciLRtidy} object is created using bpcpfit (plotstyle = "ggplot"), the plot will display using ggplot. 
If a \code{kmciLRgroup} object is created using bpcpfit (plotstyle = "standard"), the plot will display using base R.
}
\usage{

\method{plot}{kmci}(x, ...)

\method{plot}{kmciLR}(x, XLAB = "time", YLAB = "Survival", YLIM = c(0, 1), 
    ciLTY = 2, ciCOL = gray(0.8), mark.time = NULL, 
    linetype = "both", ...)

\method{lines}{kmciLR}(x, lty = c(2, 1), col = c(gray(0.8), gray(0)), 
    linetype = "both", mark.time = NULL, ...)

\method{lines}{kmci}(x, ...)

\method{plot}{kmciLRtidy}(x, ...)
                          
\method{plot}{kmciLRgroup}(x,XLAB="Time",YLAB="Survival",
 YLIM=c(0,1),ciLTY=2,ciCOL=gray(.8), linetype="both", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
\code{kmci}, \code{kmciLR}, \code{kmciLRtidy}, or \code{kmciLRgroup} object (created by functions described in \code{\link{kmtestALL}}) 
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
lines(fitBP,lwd=3,lty=1,col=gray(.5),linetype="ci")
legend(0,.2,legend=c("Kaplan-Meier","Thomas-Grunkemeier 95 pct C
I","Beta Product 95 pct CI"),
    lwd=c(1,1,3),lty=c(1,2,1),col=c(gray(0),gray(.8),gray(.5)))
data(leuk2)
## bpcpfit creates kmciLR object which is plotted with ggplot
fitBPsurv <- bpcpfit(Surv(time, status)~treatment, data=leuk2)
plot(fitBPsurv)
## this works even if no treatment variable is specified
plot(bpcpfit(Surv(time, status)~1, data=leuk2))

## if plotstyle "standard" is specified, a base R plot is produced
a <- bpcpfit(leuk2$time, leuk2$status, leuk2$treatment, plotstyle = "standard")
plot(a)

#can also create a tidy object to customize ggplot further
tidy <- tidykmciLR(fitBPsurv)
ggplot(tidy, aes(x = time, y = surv, ymin = lower, ymax = upper, col = group)) + 
  geom_line(show.legend=FALSE) + geom_ribbon(alpha = .2, aes(fill=group)) + 
  xlab("Time") + ylab("Survival") + 
  ggtitle("K-M curves with bpcp CIs") 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{survival}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

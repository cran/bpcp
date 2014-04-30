\name{bpcp2sampControl}
\alias{bpcp2sampControl}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Control function for \code{\link{bpcp2samp}}
}
\description{
Call function to change any one of options, and outputs a list with all defaults except argument that you changed. 
}
\usage{
bpcp2sampControl(Delta = 0, stype = "km", eps = 10^-8)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Delta}{
width of grouped confidence intervals, defaults to 0
}
  \item{stype}{
type of survival estimate, either "km" for Kaplan-Meier or "mue" for median unbiased estimator
}
  \item{eps}{
small value to make integration tractable
}
}
%\details{}
\value{
A list containing the 3 arguments.
}
%\references{}
%\author{}
%\note{}
\seealso{
 \code{\link{bpcp2samp}}
}
\examples{
bpcp2sampControl(Delta=1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{misc}
\name{summary.lmp}
\alias{summary.lmp}
\alias{print.summary.lmp}
\title{Summarize linear model fits with normal of order p errors}
\description{This function is the \code{summary} method for class \code{"lmp"}. This function
produces a set of results for a linear regression model.
By assuming that in a linear regression model the errors
are distributed as a Normal of order \eqn{p}{p}, we can use the function \code{\link{lmp}}.}
}
\usage{
\method{summary}{lmp}(object, ...)
\method{print}{summary.lmp}(x, ...)
}
\arguments{
\item{object}{An object of class \code{"lmp"}, a result of a call to \code{\link{lmp}}.}
\item{x}{An object of class \code{"summary.lmp"}.}
\item{...}{Further arguments passed to or from other methods.}
}
\value{
The function \code{summary} returns a list of summary statistics of the fitted linear model 
given in \code{\link{lmp}}, using the components (list elements) \code{\link{call}} and \code{\link{terms}}
from its argument, plus

\item{Call}{The matched call.}
\item{Residuals}{A summary of the vector of residuals \eqn{e_i}{e[i]}.}
\item{Coefficients}{Vector of coefficients.}
\item{Estimate of p}{An estimate of the structure parameter \eqn{p}{p}.}
\item{Power deviation of order p}{The power deviation of order \eqn{p}{p} given by
\deqn{S_p = \left[\frac{\sum e_i^p}{n-q}\right]^\frac{1}{p}}{%
Sp = (1/(n-p) Sum(e[i]^p))^1/p}
where \eqn{q}{q} is or the number of the estimated regression coefficients if \eqn{p}{p} is known,
either the number of the estimated regression coefficients plus 1 if \eqn{p}{p} is estimated.}
}
\examples{
x<-runif(30)
e<-rnormp(30,0,3,1.25)
y<-0.5+x+e
L<-lmp(y~x)
summary(L)
}
\keyword{regression}
\author{Angelo M. Mineo}

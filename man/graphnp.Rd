\name{graphnp}
\alias{graphnp}
\title{Plot of Normal of order p distributions}
\description{The function \code{graphnp} returns on the same device, marked with different colours,
from one to five normal of order \eqn{p}{p} curves.}
\usage{
graphnp(p=c(1,2,3,4,5), mu=0, sigmap=1, title="Normal of order p curves")
}
\arguments{
\item{p}{A vector of \eqn{p}{p} values. His length must be maximum five.}
\item{mu}{Value of the location parameter.}
\item{sigmap}{Value of the scale parameter.}
\item{title}{The title of the plot.}
}
\details{If one or more values of \code{p} are greater than or equal to 50, \code{graphnp}
will plot the density of an uniform distribution.
}
\value{A graphic device with till five different curves. The curves have different colours and
the device is completed by a legend.
}
\examples{
## Plot four different curves with p=1,2,3,4
## and 50 (it will plot a uniform distribution)
graphnp(c(1:4,50))
}
\keyword{aplot}
\author{Angelo M. Mineo}

\name{plot.simul.lmp}
\alias{plot.simul.lmp}
\title{Plots of the results of a simulation plan on a linear regression model}
\description{It returns the histograms of the estimates of the regression coefficients, of the scale
parameter \eqn{\sigma_p}{sigmap} and of the shape parameter \eqn{p}{p}.
}
\usage{
\method{plot}{simul.lmp}(x,...)
}
\arguments{
    \item{x}{A \code{simul.lmp} object, typically result of \code{\link{simul.lmp}}}
    \item{...}{Further arguments passed to or from other methods}
}
\value{The histograms of all the coefficients of the linear regression model and of the estimates of the scale
parameter \eqn{\sigma_p}{sigmap} and of the structure parameter \eqn{p}{p}.}
\examples{
sim<-simul.lmp(n=10,m=50,q=1,data=1.5,int=0,sigmap=1,p=3.5)
plot(sim)
}
\keyword{hplot}
\author{Angelo M. Mineo}

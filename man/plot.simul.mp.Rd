\name{plot.simul.mp}
\alias{plot.simul.mp}
\title{Plots of the results of a simulation plan on the parameters of a normal of order p distribution}
\description{It returns the histograms of the vector of means, estimates of \eqn{\mu}{Mp}, standard deviations,
estimates of \eqn{\sigma_p}{sigmap} and estimates of \eqn{p}{p}.
}
\usage{
\method{plot}{simul.mp}(x, ...)
}
\arguments{
    \item{x}{A \code{simul.mp} object, typically result of \code{\link{simul.mp}}}
    \item{...}{Further arguments passed to or from other methods}
}
\value{The histograms of the estimates of the parameters of a normal of order \eqn{p}{p} distribution.
}
\examples{
## The histograms of all the computed estimates
a<-simul.mp(100,50,mu=0,sigmap=1,p=3)
plot(a)
}
\keyword{hplot}
\author{Angelo M. Mineo}

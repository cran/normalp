\name{summary.simul.lmp}
\alias{summary.simul.lmp}
\alias{print.summary.simul.lmp}
\title{Summarize simulation results on linear regression model}
\description{This function is the \code{summary} method for class \code{"simul.lmp"}. This function
produces a set of results for a simulation plan for a linear regression model with
errors distributed as an exponential power distribution.
}
\usage{
\method{summary}{simul.lmp}(object, ...)
\method{print}{summary.simul.lmp}(x, ...)
}
\arguments{
\item{object}{An object of class \code{"simul.lmp"}, a result of a call to \code{\link{simul.lmp}}.}
\item{x}{An object of class \code{"summary.simul.lmp"}, usually a result of a call to
              \code{\link{summary.simul.lmp}}.}
\item{...}{Further arguments passed to or from other methods.}
}
\value{
This function returns this information:
\item{Results}{Table containing the simulation results.}
\item{Coefficients}{The true values of coefficients used on the simulation model.}
\item{Formula}{The used linear regression model.}
\item{Number of samples}{Number of samples generated.}
\item{Value of p}{Value of the shape parameter \eqn{p}{p} used to draw the samples.}
\item{Number of samples with problems on convergence}{If \eqn{p}{p} is estimated, we have information on
      the number of samples with problems on convergence.}
}
\examples{
ris<-simul.lmp(100,20,2,data=c(3,2),int=0,sigmap=1,p=3)
summary(ris)
}
\keyword{regression}
\author{Angelo M. Mineo}

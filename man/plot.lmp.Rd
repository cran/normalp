\name{plot.lmp}
\alias{plot.lmp}
\title{Diagnostic plots for an lmp object}
\description{This function produces four plots: a plot of residuals against fitted value, a Normal Q-Q plot,
a Normal of order \eqn{p}{p} Q-Q plot, a Scale-Location plot, with a \eqn{p}{p}-root of the standardized
residuals against the fitted values.
}
\usage{
\method{plot}{lmp}(x, ...)
}
\arguments{
    \item{x}{A \code{lmp} object, typically result of \code{\link{lmp}}.}
    \item{...}{Further arguments passed to or from other methods.}
}
\details{The standardized residuals in the Normal Q-Q plot are those of an object \code{\link{lm}};
in the Normal of order \eqn{p}{p} Q-Q plot and in the scale location plot the standardized residuals
are computed as \eqn{(e_i-m_p)/(s_p)}{(residuals-mp)/sp}.
}
\examples{
x<-1:20
z<-runif(20)
e<-rnormp(20,mu=0,sigmap=1,p=3)
y<-0.5+x+z+e
lmp.res<-lmp(y~x+z)
plot(lmp.res)
}
\keyword{hplot}
\author{Angelo M. Mineo}

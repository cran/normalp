\name{qqnormp}
\alias{qqnormp}
\alias{qqnormp.default}
\alias{qqlinep}
\title{Quantile-Quantile plot for a normal of order p distribution}
\description{The function \code{qqnormp} produces a normal of order \eqn{p}{p} Q-Q plot of the values in
\code{y}. The function \code{qqlinep} adds a line to a normal of order \eqn{p}{p} Q-Q plot going
through the first and the third quartile.
}
\usage{
qqnormp(y, ylim, p=2, main="Normal of order p Q-Q plot",
            xlab="Theoretical Quantiles", ylab="Sample Quantiles", ...)
qqlinep(y, p=2, ...)
}
\arguments{
    \item{y}{Vector of observations.}
    \item{p}{The structure parameter.}
    \item{main,xlab,ylab}{Plot labels.}
    \item{ylim, ...}{Graphical parameters}
}
\examples{
## Normal of order p Q-Q plot for a sample of 100 observations.
e<-rnormp(100,mu=0,sigmap=1,p=3)
qqnormp(e,p=3)
qqlinep(e,p=3)
}
\keyword{hplot}
\author{Angelo M. Mineo}

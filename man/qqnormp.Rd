\name{qqnormp}
\alias{qqnormp}
\alias{qqnormp.default}
\alias{qqlinep}
\title{Quantile-Quantile plot for an exponential power distribution}
\description{The function \code{qqnormp} produces an exponential power distribution Q-Q plot of the values in
\code{y}. The function \code{qqlinep} adds a line to an exponential power distribution Q-Q plot going
through the first and the third quartile.
}
\usage{
qqnormp(y, ylim, p, main, xlab, ylab, ...)
qqlinep(y, p=2, ...)
}
\arguments{
    \item{y}{Vector of observations.}
    \item{p}{The shape parameter.}
    \item{main,xlab,ylab}{Plot labels.}
    \item{ylim, ...}{Graphical parameters}
}
\examples{
## Exponential power distribution Q-Q plot for a sample of 100 observations.
e<-rnormp(100,mu=0,sigmap=1,p=3)
qqnormp(e,p=3)
qqlinep(e,p=3)
}
\keyword{hplot}
\author{Angelo M. Mineo}

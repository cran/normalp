\name{paramp}
\alias{paramp}
\alias{paramp.default}
\alias{print.paramp}
\title{Estimation of location and scale parameters}
\description{The function \code{paramp} returns a list with five elements: arithmetic mean, \eqn{M_p}{Mp},
standard deviation, \eqn{S_p}{Sp}, and shape parameter \eqn{p}{p}, estimated on a sample.
}
\usage{
paramp(x, p)
}
\arguments{
    \item{x}{A vector of observations.}
    \item{p}{If specified, the algorithm uses this value for \code{p},
                i.e. the algorithm does not use an estimate of \code{p}.}
    \item{...}{Further arguments passed to or from other methods.}
}
\value{The estimation of \eqn{\mu}{%
Mp} and \eqn{p}{p} is based on an iterative method. To show differences between the
least squares method and the maximum likelihood method, it prints out also the mean and the standard deviation.
\item{Mean}{Arithmetic mean.}
\item{Mp}{The estimated value of the location parameter.}
\item{Sd}{Standard deviation.}
\item{Sp}{The estimated value of the scale parameter.}
\item{p}{The estimated value of the shape parameter.}
\item{iter}{If its value is 1, we have had problems on convergence.}
}
\references{Mineo, A.M. (1996) \emph{La migliore combinazione delle osservazioni: curve normali di ordine p e
stimatori di norma Lp}. PhD thesis.
}
\examples{
x<-rnormp(1000,2,3,4.2)
paramp(x)
}
\keyword{univar}
\author{Angelo M. Mineo}

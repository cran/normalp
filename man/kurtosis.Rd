\name{kurtosis}
\alias{kurtosis}
\title{Indices of kurtosis}
\description{This function computes the theoretical and empirical values of three indices of kurtosis.
}
\usage{
kurtosis(x = NULL, p, value = c("estimate", "parameter"))
}
\arguments{
    \item{x}{A sample of observations.}
    \item{p}{the shape parameter.}
    \item{value}{If is set to \code{estimate}, evaluate the indices using an estimate of \eqn{p}{p}.
    Otherwise, if is set to \code{parameter} it uses the value specified in \code{p}.}
}
\value{
It returns the vector of the three indices of kurtosis \eqn{VI}{VI}, \eqn{\beta_2}{Beta2}
and \eqn{\beta_p}{Betap}.
Giving a vector as argument, it returns the estimates of the three indices, computed on the sample.
On the other hand, giving the value of the shape parameter \eqn{p}{p}, it returns the theoretical indices.
}
\references{Mineo, A.M. (1996) \emph{La migliore combinazione delle osservazioni: curve normali di ordine p e
stimatori di norma Lp}. PhD thesis.
}
\examples{
kurtosis(p=2)
x<-rnormp(50,mu=0,sigmap=2,p=1.5)
kurtosis(x,p=2)
}
\keyword{univar}
\author{Angelo M. Mineo}

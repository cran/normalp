\name{kurtosis}
\alias{kurtosis}
\title{Indexes of kurtosis}
\description{This function computes the theoretical and empirical values of three indexes of kurtosis.
}
\usage{
kurtosis(x, p=3, value="parameter")
}
\arguments{
    \item{x}{A sample of observations.}
    \item{p}{the structure parameter.}
    \item{value}{If is set to \code{estimate}, evaluate the indexes using an estimate of \eqn{p}{p}.
    Otherwise, if is set to \code{parameter} it uses the value specified in \code{p}.}
}
\value{
It returns the vector of the three indexes of kurtosis \eqn{VI}{VI}, \eqn{\beta_2}{Beta2}
and \eqn{\beta_p}{Betap}.
Giving a vector as argument, it returns the estimates of the three indexes, computed on the sample.
On the other hand, giving the value of the structure parameter \eqn{p}{p}, it returns the theoretical indexes.
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

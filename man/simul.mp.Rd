\name{simul.mp}
\alias{simul.mp}
\alias{print.simul.mp}
\alias{simul.mp.default}
\title{Simulation planning for the parameters of an exponential power distribution}
\description{This function performs a Monte Carlo simulation to compare least square estimators and
Maximum Likelihood estimators for the parameters of an exponential power distribution.
For each sample, it calls the function \code{\link{paramp}}, returning the arithmetic
means, the max-likelihood estimates of the location parameter, the standard deviations,
the max-likelihood estimates of the scale parameter and the estimates of the shape parameter.
}
\usage{
simul.mp(n, m, mu=0, sigmap=1, p=2)
}
\arguments{
    \item{n}{Sample size.}
    \item{m}{Number of samples.}
    \item{mu}{Value of the location parameter.}
    \item{sigmap}{Value of the scale parameter.}
    \item{p}{the shape parameter.}
}
\value{This function is useful to compare several kinds of estimators.
It returns an object of class \code{"simul.mp"}, a list containing the following components:
\item{dat}{A matrix \eqn{m\times 5}{m x 5} containing the results of \code{\link{paramp}} for each sample.}
\item{table}{A matrix reporting the means and the variances of the values of the five estimators.}
}
\references{Mineo, A.M. (1995) \emph{Stima dei parametri di intensit\`a e di scala di una curva normale
di ordine p (p incognito)}. Annali della Facolt\`a di Economia dell'Universit\`a di Palermo
(Area Statistico-Matematica), pp. 125-159.
}
\examples{
## Simulation plan for 100 samples of size 20, with mu=0, sigmap=1, p=3.
simul.mp(20,100,mu=0,sigmap=1,p=3)
}
\keyword{univar}
\author{Angelo M. Mineo}

\name{simul.lmp}
\alias{simul.lmp}
\title{Simulation planning for a linear regression model with errors distributed as a Normal of Order p}
\description{This function performs a Monte Carlo simulation to compare least squares estimators and
Maximum Likelihood estimators for a linear regression model with errors distributed as a normal of order
\eqn{p}{p} distribution. The regressors are drawn from a \link{Uniform} distribution.
}
\usage{
simul.lmp(n, m, q, data, int=0, sigmap=1, p=2, lp=FALSE)
}
\arguments{
    \item{n}{Sample size.}
    \item{m}{Number of samples.}
    \item{q}{Number of regressors.}
    \item{data}{A vector of coefficients.}
    \item{int}{Value of the intercept.}
    \item{sigmap}{The scale parameter.}
    \item{p}{The structure parameter.}
    \item{lp}{Logical. If \code{TRUE}, it evaluates the coefficients with p known.}
    }
\value{The function \code{simul.lmp} returns an object of class \code{"simul.lmp"}. A component of this object
is a table of means and variances of the \eqn{m}{m} estimates of the regression coefficients and
of the scale paramenter \eqn{\sigma_p}{sigmap}.
The \code{summary} shows this table and the arguments of the simulation plan. The function \code{plot}
returns the histograms of the computed estimates.
}
\references{Mineo, A.M. (1995) \emph{Stima dei parametri di regressione lineare semplice quando gli errori seguono una
distribuzione normale di ordine p (p incognito)}. Annali della Facolt\`a di Economia dell'Universit\`a
di Palermo (Area Statistico-Matematica), pp. 161-186.
}
\examples{
## Simulation of 50 samples of size 10 for a linear regression model with 1 regressor.
simul.lmp(10,50,1,data=1.5,int=1,sigmap=1,p=3,lp=FALSE)
}
\keyword{regression}
\author{Angelo M. Mineo}

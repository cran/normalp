\name{lmp}
\alias{lmp}
\alias{lmp.default}
\title{Fitted linear model with exponential power distribution errors}
\description{The function \code{lmp} is used to fit linear model. It can be used when the errors are
distributed as an exponential power distribution.
}
\usage{
lmp(formula, data, p)
}
\arguments{
\item{formula}{A symbolic description of the model to be fitted.}
\item{data}{An optional data frame containing the variables in the model.
By default the variables are taken from the environment.}
\item{p}{The shape parameter. If specified, this function estimates the parameter by using the
\eqn{L_p-norm}{Lp-norm} method.}
}
\details{To evaluate the coefficients of the linear model, \code{lmp} uses the maximum likelihood estimators.
This function can give some problems if the number of regressors is too high.}
\value{
The function \code{lmp} returns an object of \code{\link{class}} \code{"lmp"} and \code{"lm"}. The function
\code{summary} print a summary of the results.
The generic accessor functions \code{coefficients, effects, fitted.values} and
\code{residuals} extract various useful features of the value returned by \code{lmp}.
An object of \code{\link{class}} \code{"lmp"} is a list containing at least the following components:
 \item{coefficients}{A named vector of coefficients.}
 \item{residuals}{The residuals, that is responses minus fitted values.}
 \item{fitted.values}{The fitted values.}
 \item{rank}{The numeric rank of the fitted linear model.}
 \item{df.residual}{The residual degrees of freedom computed as in \code{lm}.}
 \item{call}{The matched call.}
 \item{terms}{The \code{\link{terms}} object used.}
 \item{p}{Estimate of the shape parameter computed on residuals.}
 \item{knp}{A logical parameter used by \code{summary}.}
 \item{model}{The model frame used.}
 \item{iter}{If its value is 1 we have had a difficult convergence.}
}
\references{Mineo, A.M. (1995) \emph{Stima dei parametri di regressione lineare semplice quando gli errori 
seguono una distribuzione normale di ordine p (p incognito)}. Annali della Facolt\`a di Economia dell'Universit\`a
di Palermo (Area Statistico-Matematica), pp. 161-186.
}
\examples{
e<-rnormp(n=100,mu=0,sigmap=4,p=3,method="d")
x<-runif(100)
y<-0.5+2*x+e
lmp(y~x)
}
\keyword{regression}
\author{Angelo M. Mineo}

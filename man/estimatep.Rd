\name{estimatep}
\alias{estimatep}
\title{Estimation of p}
\description{The \code{estimatep} function estimates the shape parameter \eqn{p} from a vector of observations.
}
\usage{
estimatep(x, mu, p=2, method=c("inverse","direct"))
}
\arguments{
    \item{x}{Vector of observations.}
    \item{mu}{An estimate of the location parameter.}
    \item{p}{Starting value of the shape parameter.}
    \item{method}{Method used to estimate \eqn{p}{p} from a sample.}
}
\details{The used algorithm is based on a method proposed by A.M. Mineo (1994), which uses a particular
 index of kurtosis, called \eqn{VI}{VI}

\deqn{VI=\frac{\sqrt{\Gamma(1/p)\Gamma(3/p)}}{\Gamma(2/p)}.}{%
VI=[(G(1/p) G(3/p))^{1/2}]/G(2/p).}

With \code{method} the user can choice between an inverse interpolation (faster) or a
direct solution of the equation

\deqn{\hat{VI}=\frac{\sqrt{\Gamma(1/p)\Gamma(3/p)}}{\Gamma(2/p)}.}{%
(Estimate o VI)=[(G(1/p) G(3/p))^{1/2}]/G(2/p).}
}
\value{An estimate of \eqn{p}{p} from a sample of observations.}
\references{Mineo, A.M. (1994) \emph{Un nuovo metodo di stima di p per una corretta valutazione dei
parametri di intensit\`a e di scala di una curva normale di ordine p}, 
Atti della XXXVII Riunione Scientifica della Società Italiana di Statistica, San Remo, Vol. 2, pp. 147-154.
}
\examples{
x<-rnormp(300,mu=1,sigmap=2,p=4)
p<-estimatep(x,mu=1,p=2)
p
}
\keyword{univar}
\author{Angelo M. Mineo}

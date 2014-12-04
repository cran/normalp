\name{dnormp}
\alias{dnormp}
\title{Density function of an exponential power distribution}
\description{Density function for the exponential power distribution with location parameter \code{mu},
scale parameter \code{sigmap} and shape parameter \code{p}.}
\usage{
dnormp(x, mu=0, sigmap=1, p=2, log=FALSE)
}
\arguments{
    \item{x}{Vector of quantiles.}
    \item{mu}{Vector of location parameters.}
    \item{sigmap}{Vector of scale parameters.}
    \item{p}{Shape parameter.}
    \item{log}{Logical; if TRUE, the density is given as log(density).}
}
\details{If \code{mu}, \code{sigmap} or \code{p} are not specified they assume the default values 0, 1 and 2,
respectively.
The exponential power distribution has density function

\deqn{f(x) = \frac{1}{2 p^{(1/p)} \Gamma(1+1/p) \sigma_p} e^{-\frac{|x - \mu|^p}{p \sigma_p^p}}}{%
f(x) = 1/(2 p^(1/p) Gamma(1+1/p) sigmap) exp{-|x - mu|^p/(p sigmap^p)}}

where \eqn{\mu}{mu} is the location parameter, \eqn{\sigma_p}{sigmap} the scale parameter and \eqn{p}{p} the
shape parameter.
When \eqn{p=2}{p=2} the exponential power distribution becomes the Normal Distribution, when
\eqn{p=1}{p=1} the exponential power distribution becomes the Laplace Distribution, when
\eqn{p\rightarrow\infty}{p->infinity} the exponential power distribution becomes the Uniform Distribution.
}
\value{\code{dnormp} gives the density function of an exponential power distribution.}
\seealso{\code{\link{Normal}} for the Normal distribution, \code{\link{Uniform}} for the Uniform distribution,
and \code{\link{Special}} for the Gamma function. }
\examples{
## Compute the density for a vector x with mu=0, sigmap=1 and p=1.5
## At the end we have the graph of the exponential power distribution 
## density function with p=1.5
x <- c(-1, 1)
f <- dnormp(x, p=1.5)
print(f)
plot(function(x) dnormp(x, p=1.5) , -4, 4,
          main = "Exponential power distribution density function (p=1.5)", ylab="f(x)")
}
\keyword{distribution}
\author{Angelo M. Mineo}

\name{rnormp}
\alias{rnormp}
\title{Pseudo-random numbers from an exponential power distribution}
\description{Generation of pseudo-random numbers from an exponential power distribution
with location parameter \code{mu}, scale parameter \code{sigmap} and shape parameter \code{p}.
}
\usage{rnormp(n, mu = 0, sigmap = 1, p = 2, method = c("def", "chiodi"))}
\arguments{
    \item{n}{Number of observations.}
    \item{mu}{Vector of location parameters.}
    \item{sigmap}{Vector of scale parameters.}
    \item{p}{Shape parameter.}
    \item{method}{If is set to the default method "\code{def}", it uses the method based on the
transformation of a Gamma random variable. If set to "\code{chiodi}", it uses an algorithm based on
a generalization of the Marsaglia formula to generate pseudo-random numbers from a normal distribution.
The default method "\code{def}" is faster than the "\code{chiodi}" one
(this one is introduced only for "historical" purposes).}
}
\details{If \code{mu}, \code{sigmap} or \code{p} are not specified they assume the default values 0, 1 and 2,
respectively.
The exponential power distribution has density function

\deqn{f(x) = \frac{1}{2 p^{(1/p)} \Gamma(1+1/p) \sigma_p} e^{- \frac{|x - \mu|^p}{p \sigma_p^p}}}{%
f(x) = 1/(2 p^(1/p) Gamma(1+1/p) sigmap) exp{-|x - mu|^p/(p sigmap^p)}}

where \eqn{\mu}{mu} is the location parameter, \eqn{\sigma_p}{sigmap} the scale parameter and \eqn{p}{p} the
shape parameter.
When \eqn{p=2}{p=2} the exponential power distribution becomes the Normal Distribution, when
\eqn{p=1}{p=1} the exponential power distribution becomes the Laplace Distribution, when
\eqn{p\rightarrow\infty}{p->infinity} the exponential power distribution becomes the Uniform Distribution.
}
\value{\code{rnormp} gives a vector of \eqn{n}{n} pseudo-random numbers from an exponential power distribution.}
\seealso{\code{\link{Normal}} for the Normal distribution, \code{\link{Uniform}} for the Uniform distribution,
\code{\link{Special}} for the Gamma function and \code{\link{.Random.seed}} for the random number generation.}
\references{
Chiodi, M. (1986) \emph{Procedures for generating pseudo-random numbers from a normal distribution of order p (p>1)},
Statistica Applicata, 1, pp. 7-26.

Marsaglia, G. and Bray, T.A. (1964) \emph{A convenient method for generating normal variables},
SIAM rev., 6, pp. 260-264.
}
\examples{
## Generate a random sample x from an exponential power distribution
## At the end we have the histogram of x
x <- rnormp(1000, 1, 2, 1.5)
hist(x, main="Histogram of the random sample")
}
\keyword{distribution}
\author{Angelo M. Mineo}

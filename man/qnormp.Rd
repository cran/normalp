\name{qnormp}
\alias{qnormp}
\title{Quantiles of an exponential power distribution}
\description{Quantiles for the exponential power distribution with location parameter \code{mu},
scale parameter \code{sigmap} and shape parameter \code{p}.}
\usage{
qnormp(pr, mu=0, sigmap=1, p=2, lower.tail=TRUE, log.pr=FALSE)
}
\arguments{
    \item{pr}{Vector of probabilities.}
    \item{mu}{Vector of location parameters.}
    \item{sigmap}{Vector of scale parameters.}
    \item{p}{Shape parameter.}
    \item{lower.tail}{Logical; if TRUE (default), probabilities are \eqn{P [X\leq x]}{%
    P[X <= x]}, otherwise, \eqn{P[X>x]}{P[X>x]}.}
    \item{log.pr}{Logical; if TRUE, probabilities \eqn{pr}{pr} are given as \eqn{log(pr)}{log(pr)}.}
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
\value{\code{qnormp} gives the quantiles of an exponential power distribution.}
\seealso{\code{\link{Normal}} for the Normal distribution, \code{\link{Uniform}} for the Uniform distribution, and \code{\link{Special}} for the Gamma function. }
\examples{
## Compute the quantiles for a vector of probabilities x
## with mu=1, sigmap=2 and p=1.5
x <- 0.3
q <- qnormp(x, 1, 2, 1.5)
q
}
\keyword{distribution}
\author{Angelo M. Mineo}

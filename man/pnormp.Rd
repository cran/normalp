\name{pnormp}
\alias{pnormp}
\title{Probability function of a Normal of Order p Distribution}
\description{Probability function for the normal of order \eqn{p}{p} distribution with location parameter
\code{mu}, scale parameter \code{sigmap} and structure parameter \code{p}.
}
\usage{
pnormp(q, mu=0, sigmap=1, p=2, lower.tail=TRUE, log.pr=FALSE)
}
\arguments{
    \item{q}{Vector of quantiles.}
    \item{mu}{Vector of location parameters.}
    \item{sigmap}{Vector of scale parameters.}
    \item{p}{Structure parameter.}
    \item{lower.tail}{Logical; if TRUE (default), probabilities are \eqn{P [X\leq x]}{%
    P[X <= x]}, otherwise, \eqn{P[X>x]}{P[X>x]}.}
    \item{log.pr}{Logical; if TRUE, probabilities \eqn{pr}{pr} are given as \eqn{log(pr)}{log(pr)}.}
}
\details{If \code{mu}, \code{sigmap} or \code{p} are not specified they assume the default values 0, 1 and 2,
respectively.
The normal of order \eqn{p}{p} distribution has density function

\deqn{f(x) = \frac{1}{2 p^{(1/p)} \Gamma(1+1/p) \sigma_p} e^{-\frac{|x - \mu|^p}{p \sigma_p^p}}}{%
f(x) = 1/(2 p^(1/p) Gamma(1+1/p) sigmap) exp{-|x - mu|^p/(p sigmap^p)}}

where \eqn{\mu}{mu} is the location parameter, \eqn{\sigma_p}{sigmap} the scale parameter and \eqn{p}{p} the
structure parameter.
When \eqn{p=2}{p=2} the Normal of Order \eqn{p}{p} Distribution becomes the Normal (Gaussian) Distribution,
when \eqn{p=1}{p=1} the Normal of Order \eqn{p}{p} Distribution becomes the Laplace Distribution, when
\eqn{p\rightarrow\infty}{p->infinity} the Normal of Order \eqn{p}{p} Distribution becomes the Uniform Distribution.
}
\value{\code{pnormp} gives the probability of a normal of order \eqn{p}{p} distribution.}
\seealso{\code{\link{Normal}} for the Normal distribution, \code{\link{Uniform}} for the Uniform distribution, and \code{\link{Special}} for the Gamma function. }
\examples{
## Compute the distribution function for a vector x with mu=0, sigmap=1 and p=1.5
## At the end we have the graph of the Normal of order 1.5 distribution function.
x <- c(-1, 1)
pr <- pnormp(x, p=1.5)
print(pr)
plot(function(x) pnormp(x, p=1.5), -4, 4,
          main = "Normal of order p Distribution Function (p=1.5)", ylab="F(x)")
}
\keyword{distribution}
\author{Angelo M. Mineo}

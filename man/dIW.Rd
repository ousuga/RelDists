% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dIW.R
\name{dIW}
\alias{dIW}
\alias{pIW}
\alias{qIW}
\alias{rIW}
\alias{hIW}
\title{The Inverse Weibull distribution}
\usage{
dIW(x, mu, sigma, log = FALSE)

pIW(q, mu, sigma, lower.tail = TRUE, log.p = FALSE)

qIW(p, mu, sigma, lower.tail = TRUE, log.p = FALSE)

rIW(n, mu, sigma)

hIW(x, mu, sigma)
}
\arguments{
\item{x, q}{vector of quantiles.}

\item{mu}{scale parameter.}

\item{sigma}{shape parameters.}

\item{log, log.p}{logical; if TRUE, probabilities p are given as log(p).}

\item{lower.tail}{logical; if TRUE (default), probabilities are P[X <= x], otherwise, P[X > x].}

\item{p}{vector of probabilities.}

\item{n}{number of observations.}
}
\value{
\code{dIW} gives the density, \code{pIW} gives the distribution 
function, \code{qIW} gives the quantile function, \code{rIW}
generates random deviates and \code{hIW} gives the hazard function.
}
\description{
Density, distribution function, quantile function, 
random generation and hazard function for the inverse weibull distribution with
parameters \code{mu} and \code{sigma}.
}
\details{
The inverse weibull distribution with parameters \code{mu} and
\code{sigma} has density given by

\eqn{f(x) = \mu \sigma x^{-\sigma-1} \exp(-\mu x^{-\sigma})}

for \eqn{x > 0}, \eqn{\mu > 0} and \eqn{\sigma > 0}
}
\examples{
# The probability density function
curve(dIW(x, mu=1, sigma=2), from=0, to=10,
      col="red", las=1, ylab="f(x)")

# The cumulative distribution and the Reliability function
par(mfrow=c(1, 2))
curve(pIW(x, mu=1, sigma=2),
      from=0, to=10,  col="red", las=1, ylab="F(x)")
curve(pIW(x, mu=1, sigma=2, lower.tail=FALSE),
      from=0, to=10, col="red", las=1, ylab="R(x)")
            
# The quantile function
p <- seq(from=0, to=0.99, length.out=100)
plot(x=qIW(p, mu=1, sigma=2), y=p, xlab="Quantile",
  las=1, ylab="Probability")
curve(pIW(x, mu=1, sigma=2), from=0, add=TRUE, col="red")
  
# The random function
hist(rIW(n=1000, mu=1, sigma=2), freq=FALSE, xlim=c(0, 40),
  xlab="x", las=1, main="")
curve(dIW(x, mu=1, sigma=2), from=0, add=TRUE, col="red")

# The Hazard function
par(mfrow=c(1, 1))
curve(hIW(x, mu=1, sigma=2), from=0, to=15,
   col="red", ylab="Hazard function", las=1)

}
\references{
\insertRef{almalki2014modifications}{RelDists}

\insertRef{drapella1993complementary}{RelDists}
}
\author{
Freddy Hernandez, \email{fhernanb@unal.edu.co}
}

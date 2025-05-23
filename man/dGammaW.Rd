% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dGammaW.R
\name{dGammaW}
\alias{dGammaW}
\alias{pGammaW}
\alias{qGammaW}
\alias{rGammaW}
\alias{hGammaW}
\title{The Gamma Weibull distribution}
\usage{
dGammaW(x, mu, sigma, nu, log = FALSE)

pGammaW(q, mu, sigma, nu, lower.tail = TRUE, log.p = FALSE)

qGammaW(p, mu, sigma, nu, lower.tail = TRUE, log.p = FALSE)

rGammaW(n, mu, sigma, nu)

hGammaW(x, mu, sigma, nu)
}
\arguments{
\item{x, q}{vector of quantiles.}

\item{mu}{parameter.}

\item{sigma}{parameter.}

\item{nu}{parameter.}

\item{log, log.p}{logical; if TRUE, probabilities p are given as log(p).}

\item{lower.tail}{logical; if TRUE (default), probabilities are 
P[X <= x], otherwise, P[X > x].}

\item{p}{vector of probabilities.}

\item{n}{number of observations.}
}
\value{
\code{dGammaW} gives the density, \code{pGammaW} gives the distribution 
function, \code{qGammaW} gives the quantile function, \code{rGammaW}
generates random deviates and \code{hGammaW} gives the hazard function.
}
\description{
Density, distribution function, quantile function, 
random generation and hazard function for the Gamma Weibull distribution 
with parameters \code{mu}, \code{sigma}, \code{nu} and \code{tau}.
}
\details{
The Gamma Weibull Distribution with parameters \code{mu}, 
\code{sigma} and \code{nu} has density given by

\eqn{f(x)= \frac{\sigma \mu^{\nu}}{\Gamma(\nu)} x^{\nu \sigma - 1} \exp(-\mu x^\sigma),}

for \eqn{x > 0}, \eqn{\mu > 0}, \eqn{\sigma \geq 0} and \eqn{\nu > 0}.
}
\examples{
old_par <- par(mfrow = c(1, 1)) # save previous graphical parameters

## The probability density function 
curve(dGammaW(x, mu = 0.5, sigma = 2, nu=1), from = 0, to = 6, 
      col = "red", las = 1, ylab = "f(x)")

## The cumulative distribution and the Reliability function
par(mfrow = c(1, 2))
curve(pGammaW(x, mu = 0.5, sigma = 2, nu=1), from = 0, to = 3, 
ylim = c(0, 1), col = "red", las = 1, ylab = "F(x)")
curve(pGammaW(x, mu = 0.5, sigma = 2, nu=1, lower.tail = FALSE), 
from = 0, to = 3, ylim = c(0, 1), col = "red", las = 1, ylab = "R(x)")

## The quantile function
p <- seq(from = 0, to = 0.99999, length.out = 100)
plot(x = qGammaW(p = p, mu = 0.5, sigma = 2, nu=1), y = p, 
xlab = "Quantile", las = 1, ylab = "Probability")
curve(pGammaW(x, mu = 0.5, sigma = 2, nu=1), from = 0, add = TRUE, 
col = "red")

## The random function
hist(rGammaW(1000, mu = 0.5, sigma = 2, nu=1), freq = FALSE, xlab = "x", 
ylim = c(0, 1), las = 1, main = "")
curve(dGammaW(x, mu = 0.5, sigma = 2, nu=1),  from = 0, add = TRUE, 
col = "red", ylim = c(0, 1))

## The Hazard function
par(mfrow=c(1,1))
curve(hGammaW(x, mu = 0.5, sigma = 2, nu=1), from = 0, to = 2, 
ylim = c(0, 1), col = "red", ylab = "Hazard function", las = 1)

par(old_par) # restore previous graphical parameters
}
\references{
Almalki, S. J., & Nadarajah, S. (2014). Modifications of the 
Weibull distribution: A review. Reliability Engineering & 
System Safety, 124, 32-55.

Stacy, E. W. (1962). A generalization of the gamma 
distribution. The Annals of mathematical statistics, 1187-1192.
}
\seealso{
\link{GammaW}
}
\author{
Johan David Marin Benjumea, \email{johand.marin@udea.edu.co}
}

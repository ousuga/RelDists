% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dKumIW.R
\name{dKumIW}
\alias{dKumIW}
\alias{pKumIW}
\alias{qKumIW}
\alias{rKumIW}
\alias{hKumIW}
\title{The Kumaraswamy Inverse Weibull distribution}
\usage{
dKumIW(x, mu, sigma, nu, log = FALSE)

pKumIW(q, mu, sigma, nu, lower.tail = TRUE, log.p = FALSE)

qKumIW(p, mu, sigma, nu, lower.tail = TRUE, log.p = FALSE)

rKumIW(n, mu, sigma, nu)

hKumIW(x, mu, sigma, nu)
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
\code{dKumIW} gives the density, \code{pKumIW} gives the distribution 
function, \code{qKumIW} gives the quantile function, \code{rKumIW}
generates random deviates and \code{hKumIW} gives the hazard function.
}
\description{
Density, distribution function, quantile function, 
random generation and hazard function for the Kumaraswamy Inverse Weibull distribution 
with parameters \code{mu}, \code{sigma} and \code{nu}.
}
\details{
The Kumaraswamy Inverse Weibull Distribution with parameters \code{mu}, 
\code{sigma} and \code{nu} has density given by

\eqn{f(x)= \mu \sigma \nu x^{-\mu - 1} \exp{- \sigma x^{-\mu}} (1 - \exp{- \sigma x^{-\mu}})^{\nu - 1},}

for \eqn{x > 0}, \eqn{\mu > 0}, \eqn{\sigma > 0} and \eqn{\nu > 0}.
}
\examples{
old_par <- par(mfrow = c(1, 1)) # save previous graphical parameters

## The probability density function 
par(mfrow = c(1, 1))
curve(dKumIW(x, mu = 1.5, sigma=  1.5, nu = 1), from = 0, to = 8.5, 
      col = "red", las = 1, ylab = "f(x)")

## The cumulative distribution and the Reliability function
par(mfrow = c(1, 2))
curve(pKumIW(x, mu = 1.5, sigma=  1.5, nu = 1), from = 0, to = 8.5, 
      ylim = c(0, 1), col = "red", las = 1, ylab = "F(x)")
curve(pKumIW(x, mu = 1.5, sigma=  1.5, nu = 1, lower.tail = FALSE), 
      from = 0, to = 6, ylim = c(0, 1), col = "red", las = 1, ylab = "R(x)")

## The quantile function
p <- seq(from = 0, to = 0.99999, length.out = 100)
plot(x = qKumIW(p=p, mu = 1.5, sigma=  1.5, nu = 10), y = p, 
     xlab = "Quantile", las = 1, ylab = "Probability")
curve(pKumIW(x, mu = 1.5, sigma=  1.5, nu = 10), from = 0, add = TRUE, 
      col = "red")

## The random function
hist(rKumIW(1000, mu = 1.5, sigma=  1.5, nu = 5), freq = FALSE, xlab = "x", 
     las = 1, ylim = c(0, 1.5), main = "")
curve(dKumIW(x, mu = 1.5, sigma=  1.5, nu = 5), from = 0, to =8, add = TRUE, 
      col = "red")

## The Hazard function
par(mfrow=c(1,1))
curve(hKumIW(x, mu = 1.5, sigma=  1.5, nu = 1), from = 0, to = 3, 
      ylim = c(0, 0.7), col = "red", ylab = "Hazard function", las = 1)

par(old_par) # restore previous graphical parameters
}
\references{
Almalki, S. J., & Nadarajah, S. (2014). Modifications of the 
Weibull distribution: A review. Reliability Engineering & 
System Safety, 124, 32-55.

Shahbaz, M. Q., Shahbaz, S., & Butt, N. S. (2012). 
The Kumaraswamy Inverse Weibull Distribution. 
Pakistan journal of statistics and operation research, 479-489.
}
\seealso{
\link{KumIW}
}
\author{
Johan David Marin Benjumea, \email{johand.marin@udea.edu.co}
}

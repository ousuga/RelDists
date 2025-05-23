% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/QXGP.R
\name{QXGP}
\alias{QXGP}
\title{The Quasi XGamma Poisson family}
\usage{
QXGP(mu.link = "log", sigma.link = "log", nu.link = "log")
}
\arguments{
\item{mu.link}{defines the mu.link, with "log" link as the default for the mu parameter.}

\item{sigma.link}{defines the sigma.link, with "log" link as the default for the sigma.}

\item{nu.link}{defines the nu.link, with "log" link as the default for the nu parameter.}
}
\value{
Returns a gamlss.family object which can be used to fit a QXGP distribution in the \code{gamlss()} function.
}
\description{
The Quasi XGamma Poisson family
}
\details{
The Quasi XGamma Poisson distribution with parameters \code{mu}, 
\code{sigma} and \code{nu} has density given by

\eqn{f(x)= K(\mu, \sigma, \nu)(\frac {\sigma^{2} x^{2}}{2} + \mu)
 exp(\frac{\nu exp(-\sigma x)(1 + \mu + \sigma x + \frac {\sigma^{2}x^{2}}{2})}{1+\mu} - \sigma x),}

for \eqn{x > 0}, \eqn{\mu> 0}, \eqn{\sigma> 0}, \eqn{\nu> 1}.

where

\eqn{K(\mu, \sigma, \nu) = \frac{\nu \sigma}{(exp(\nu)-1)(1+\mu)}}
}
\examples{
# Example 1
# Generating some random values with
# known mu, sigma and nu
y <- rQXGP(n=200, mu=4, sigma=2, nu=3)

# Fitting the model
require(gamlss)

mod <- gamlss(y~1, sigma.fo=~1, nu.fo=~1, family='QXGP',
              control=gamlss.control(n.cyc=5000, trace=FALSE))

# Extracting the fitted values for mu, sigma and nu
# using the inverse link function
exp(coef(mod, what='mu'))
exp(coef(mod, what='sigma'))
exp(coef(mod, what='nu'))

# Example 2
# Generating random values under some model
n <- 2000
x1 <- runif(n, min=0.4, max=0.6)
x2 <- runif(n, min=0.4, max=0.6)
mu <- exp(-2.19 + 3 * x1)
sigma <- exp(1 - 2 * x2)
nu <- 1
x <- rQXGP(n=n, mu, sigma, nu)

mod <- gamlss(x~x1, sigma.fo=~x2, nu.fo=~1, family=QXGP,
              control=gamlss.control(n.cyc=5000, trace=FALSE))

coef(mod, what="mu")
coef(mod, what="sigma")
exp(coef(mod, what="nu"))
}
\references{
Sen, S., Korkmaz, M. Ç., & Yousof, H. M. (2018). 
The quasi XGamma-Poisson distribution: properties and
application. Istatistik Journal of The Turkish Statistical 
Association, 11(3), 65-76.
}
\seealso{
\link{dQXGP}
}
\author{
Amylkar Urrea Montoya, \email{amylkar.urrea@udea.edu.co}
}

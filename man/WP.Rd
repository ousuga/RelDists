% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/WP.R
\name{WP}
\alias{WP}
\title{The Weibull Poisson family}
\usage{
WP(mu.link = "log", sigma.link = "log", nu.link = "log")
}
\arguments{
\item{mu.link}{defines the mu.link, with "log" link as the default for the mu parameter.}

\item{sigma.link}{defines the sigma.link, with "log" link as the default for the sigma.}

\item{nu.link}{defines the nu.link, with "log" link as the default for the nu parameter.}
}
\value{
Returns a gamlss.family object which can be used to fit a WP distribution in the \code{gamlss()} function.
}
\description{
The Weibull Poisson family
}
\details{
The Weibull Poisson distribution with parameters \code{mu}, 
\code{sigma} and \code{nu} has density given by

\eqn{f(x) = \frac{\mu \sigma \nu e^{-\nu}} {1-e^{-\nu}} x^{\mu-1} exp({-\sigma x^{\mu}+\nu exp({-\sigma} x^{\mu}) }),}

for x > 0.
}
\examples{
# Example 1
# Generating some random values with
# known mu, sigma and nu
y <- rWP(n=300, mu=1.5, sigma=0.5, nu=0.5)

# Fitting the model
require(gamlss)

mod <- gamlss(y~1, sigma.fo=~1, nu.fo=~1, family='WP',
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
mu <- exp(-1.3 + 3 * x1)
sigma <- exp(0.69 - 2 * x2)
nu <- 0.5
x <- rWP(n=n, mu, sigma, nu)

mod <- gamlss(x~x1, sigma.fo=~x2, nu.fo=~1, family=WP,
              control=gamlss.control(n.cyc=5000, trace=FALSE))

coef(mod, what="mu")
coef(mod, what="sigma")
exp(coef(mod, what="nu"))
}
\references{
\insertRef{almalki2014modifications}{RelDists}

\insertRef{Wanbo2012}{RelDists}
}
\seealso{
\link{dWP}
}
\author{
Amylkar Urrea Montoya, \email{amylkar.urrea@udea.edu.co}
}

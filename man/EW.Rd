% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/EW.R
\name{EW}
\alias{EW}
\title{The Exponentiated Weibull family}
\usage{
EW(mu.link = "log", sigma.link = "log", nu.link = "log")
}
\arguments{
\item{mu.link}{defines the mu.link, with "log" link as the default for the mu parameter.}

\item{sigma.link}{defines the sigma.link, with "log" link as the default for the sigma.}

\item{nu.link}{defines the nu.link, with "log" link as the default for the nu parameter.}
}
\value{
Returns a gamlss.family object which can be used to fit a EW distribution in the \code{gamlss()} function.
}
\description{
The Exponentiated Weibull distribution
}
\details{
The Exponentiated Weibull Distribution with parameters \code{mu}, 
\code{sigma} and \code{nu} has density given by

\eqn{f(x)=\nu \mu \sigma x^{\sigma-1} \exp(-\mu x^\sigma) (1-\exp(-\mu x^\sigma))^{\nu-1},}

for x > 0.
}
\examples{
# Example 1
# Generating some random values with
# known mu, sigma and nu
# Will not be run this example because high number is cycles
# is needed in order to get good estimates
\dontrun{
y <- rEW(n=100, mu=2, sigma=1.5, nu=0.5)

# Fitting the model
require(gamlss)
mod <- gamlss(y~1, sigma.fo=~1, nu.fo=~1, family='EW',
              control=gamlss.control(n.cyc=5000, trace=FALSE))

# Extracting the fitted values for mu, sigma and nu
# using the inverse link function
exp(coef(mod, what='mu'))
exp(coef(mod, what='sigma'))
exp(coef(mod, what='nu'))
}

# Example 2
# Generating random values under some model
# Will not be run this example because high number is cycles
# is needed in order to get good estimates
\dontrun{
n <- 200
x1 <- rpois(n, lambda=2)
x2 <- runif(n)
mu <- exp(2 + -3 * x1)
sigma <- exp(3 - 2 * x2)
nu <- 2
x <- rEW(n=n, mu, sigma, nu)

mod <- gamlss(x~x1, sigma.fo=~x2, nu.fo=~1, family=EW, 
              control=gamlss.control(n.cyc=5000, trace=FALSE))

coef(mod, what="mu")
coef(mod, what="sigma")
exp(coef(mod, what="nu"))
}
}
\seealso{
\link{dEW}
}

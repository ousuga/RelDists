# Example 1
# Generating some random values with
# known mu, sigma, nu and tau
\donttest{{
y <- rFPEGG(n=200, mu=0.1, sigma=0.8, nu=10, tau=1.5)

# Fitting the model
require(gamlss)

mod <- gamlss(y~1, sigma.fo=~1, nu.fo=~1, tau.fo=~1, family='FPEGG',
              control=gamlss.control(n.cyc=300, trace=FALSE))

# Extracting the fitted values for mu, sigma, nu and tau
# using the inverse link function
exp(coef(mod, what='mu'))
exp(coef(mod, what='sigma'))
exp(coef(mod, what='nu'))
exp(coef(mod, what='tau'))
}

# Example 2
# Generating random values under some model
\donttest{
n <- 200
x1 <- runif(n, min=0.4, max=0.6)
x2 <- runif(n, min=0.4, max=0.6)
mu <- exp(-0.8 + -3 * x1)
sigma <- exp(0.77 - 2 * x2)
nu <- 10
tau <- 1.5
x <- rFPEGG(n=n, mu, sigma, nu, tau)

mod <- gamlss(x~x1, sigma.fo=~x2, nu.fo=~1, tau.fo=~1, family=FPEGG,
              control=gamlss.control(n.cyc=300, trace=FALSE))

coef(mod, what="mu")
coef(mod, what="sigma")
exp(coef(mod, what="nu"))
exp(coef(mod, what="tau"))
}
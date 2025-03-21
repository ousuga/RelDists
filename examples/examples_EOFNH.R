# Example 1
# Generating some random values with
# known mu, sigma, nu and tau
set.seed(123)
y <- rEOFNH(n=100, mu=1, sigma=2.1, nu=0.8, tau=1)

# Fitting the model
require(gamlss)

mod <- gamlss(y~1, sigma.fo=~1, nu.fo=~1, tau.fo=~1, family=EOFNH,
              control=gamlss.control(n.cyc=5000, trace=FALSE))

# Extracting the fitted values for mu, sigma, nu and tau
# using the inverse link function
exp(coef(mod, what="mu"))
exp(coef(mod, what="sigma"))
exp(coef(mod, what="nu"))
exp(coef(mod, what="tau"))

# Example 2
# Generating random values under the model
n <- 100
x1 <- runif(n)
x2 <- runif(n)
mu <- exp(0.5 - 1.2 * x1)
sigma <- 2.1
nu <- 0.8
tau <- 1
y <- rEOFNH(n=n, mu, sigma, nu, tau)

mod <- gamlss(y~x1, sigma.fo=~1, nu.fo=~1, tau.fo=~1, family=EOFNH,
              control=gamlss.control(n.cyc=5000, trace=FALSE))

coef(mod, what="mu")
exp(coef(mod, what="sigma"))
exp(coef(mod, what="nu"))
exp(coef(mod, what="tau"))

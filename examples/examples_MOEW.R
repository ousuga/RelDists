# Example 1
# Generating some random values with
# known mu, sigma and nu
set.seed(12345)
y <- rMOEW(n=100, mu=20, sigma=2, nu=1)

# Fitting the model
require(gamlss)

mod <- gamlss(y~1, sigma.fo=~1, nu.fo=~1, family=MOEW,
              control=gamlss.control(n.cyc=5000, trace=TRUE))

# Extracting the fitted values for mu, sigma and nu
# using the inverse link function
exp(coef(mod, what="mu"))
exp(coef(mod, what="sigma"))
exp(coef(mod, what="nu"))

# Example 2
# Generating random values under some model
n <- 500
x1 <- runif(n, min=0.2, max=0.6)
x2 <- runif(n, min=0.2, max=0.6)
mu <- exp(-1.20 + 3 * x1)
sigma <- exp(0.84 - 2 * x2)
nu <- 1

set.seed(1234)
y <- rMOEW(n=n, mu, sigma, nu)

mod <- gamlss(y~x1, sigma.fo=~x2, nu.fo=~1, family=MOEW,
              control=gamlss.control(n.cyc=5000, trace=FALSE))

coef(mod, what="mu")
coef(mod, what="sigma")
exp(coef(mod, what="nu"))

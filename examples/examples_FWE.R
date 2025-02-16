# Example 1
# Generating some random values with
# known mu and sigma
y <- rFWE(n=100, mu=0.75, sigma=1.3)

# Fitting the model
require(gamlss)
mod <- gamlss(y~1, sigma.fo=~1, family="FWE")

# Extracting the fitted values for mu and sigma
# using the inverse link function
exp(coef(mod, what="mu"))
exp(coef(mod, what="sigma"))

# Example 2
# Generating random values under some model
n <- 200
x1 <- runif(n)
x2 <- runif(n)
mu <- exp(1.21 - 3 * x1)
sigma <- exp(1.26 - 2 * x2)
y <- rFWE(n=n, mu=mu, sigma=sigma)

mod <- gamlss(y~x1, sigma.fo=~x2, family=FWE)

coef(mod, what="mu")
coef(mod, what="sigma")

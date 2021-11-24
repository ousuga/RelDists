# Example 1
# Generating some random values with
# known mu, sigma and nu
y <- rLIN(n=200, mu=2)

# Fitting the model
require(gamlss)
mod <- gamlss(y ~ 1, family="LIN")

# Extracting the fitted values for mu
# using the inverse link function
exp(coef(mod, what='mu'))

# Example 2
# Generating random values under some model
n <- 100
x1 <- runif(n=n)
x2 <- runif(n=n)
eta <- 1 + 3 * x1 - 2 * x2
mu <- exp(eta)
y <- rLIN(n=n, mu=mu)

mod <- gamlss(y ~ x1 + x2, family=LIN)

coef(mod, what='mu')


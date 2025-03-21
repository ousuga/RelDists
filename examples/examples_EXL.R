# Example 1
# Generating some random values with
# known mu and sigma
y <- rEXL(n=1000, mu=0.75, sigma=1.3)

# Fitting the model
require(gamlss)
mod <- gamlss(y~1, sigma.fo=~1, family="EXL")

# Extracting the fitted values for mu and sigma
# using the inverse link function
exp(coef(mod, what="mu"))
exp(coef(mod, what="sigma"))


# Example 2
# Generating random values under some model
n <- 1000
x1 <- runif(n)
x2 <- runif(n)
mu <- exp(1.45 - 3 * x1)
sigma <- exp(2 - 1.5 * x2)
y <- rEXL(n=n, mu=mu, sigma=sigma)

mod <- gamlss(y~x1, sigma.fo=~x2, family=EXL)

summary(mod)

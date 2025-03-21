# Example 1
# Generating some random values with
# known mu, sigma and nu
y <- rWP(n=3000, mu=1.5, sigma=0.5, nu=0.5)

# Fitting the model
require(gamlss)

mod1 <- gamlss(y~1, sigma.fo=~1, nu.fo=~1, family=WP,
              control=gamlss.control(n.cyc=5000, trace=FALSE))

# Extracting the fitted values for mu, sigma and nu
# using the inverse link function
exp(coef(mod1, what="mu"))
exp(coef(mod1, what="sigma"))
exp(coef(mod1, what="nu"))

# Example 2
# Generating random values under some model

# A function to simulate a data set with Y ~ WP
gendat <- function(n) {
  x1 <- runif(n)
  x2 <- runif(n)
  mu <- exp(-1.3 + 3.1 * x1)
  sigma <- exp(0.9 - 3.2 * x2)
  nu <- 0.5
  y <- rWP(n=n, mu, sigma, nu)
  data.frame(y=y, x1=x1, x2)
}

set.seed(1234)
dat <- gendat(n=100)

# Fitting the model
mod2 <- NULL
mod2 <- gamlss(y~x1, sigma.fo=~x2, nu.fo=~1, 
               family=WP, data=dat,
              control=gamlss.control(n.cyc=5000, trace=FALSE))

coef(mod2, what="mu")
coef(mod2, what="sigma")
exp(coef(mod2, what="nu"))


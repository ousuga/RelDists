# Example 1
# Generating some random values with
# known mu, sigma, nu and tau
\donttest{
set.seed(123456)
y <- rEGG(n=100, mu=0.1, sigma=0.8, nu=10, tau=1.5)

# Fitting the model
require(gamlss)

mod <- gamlss(y~1, sigma.fo=~1, nu.fo=~1, tau.fo=~1, 
              family=EGG,
              control=gamlss.control(n.cyc=500, trace=FALSE))

# Extracting the fitted values for mu, sigma, nu and tau
# using the inverse link function
exp(coef(mod, what="mu"))
exp(coef(mod, what="sigma"))
exp(coef(mod, what="nu"))
exp(coef(mod, what="tau"))
}

# Example 2
# Generating random values under some model

\donttest{
# A function to simulate a data set with Y ~ EGG
gendat <- function(n) {
  x1 <- runif(n)
  x2 <- runif(n)
  mu    <- exp(-0.8 + -3 * x1)
  sigma <- exp(0.77 - 2 * x2)
  nu    <- 10
  tau   <- 1.5
  y <- rEGG(n=n, mu=mu, sigma=sigma, nu=nu, tau)
  data.frame(y=y, x1=x1, x2=x2)
}

set.seed(12345)
dat <- gendat(n=200)

mod <- gamlss(y~x1, sigma.fo=~x2, nu.fo=~1, tau.fo=~1, 
              family=EGG, data=dat,
              control=gamlss.control(n.cyc=500, trace=FALSE))

coef(mod, what="mu")
coef(mod, what="sigma")
exp(coef(mod, what="nu"))
exp(coef(mod, what="tau"))
}
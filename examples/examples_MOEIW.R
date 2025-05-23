# Example 1
# Generating some random values with
# known mu, sigma and nu
set.seed(123456)
y <- rMOEIW(n=100, mu=0.6, sigma=1.7, nu=0.3)

# Fitting the model
require(gamlss)

mod <- gamlss(y~1, sigma.fo=~1, nu.fo=~1, family="MOEIW",
              control=gamlss.control(n.cyc=5000, trace=FALSE))

# Extracting the fitted values for mu, sigma and nu
# using the inverse link function
exp(coef(mod, what="mu"))
exp(coef(mod, what="sigma"))
exp(coef(mod, what="nu"))

# Example 2
# Generating random values under some model

# A function to simulate a data set with Y ~ MOEIW
gendat <- function(n) {
  x1 <- runif(n)
  x2 <- runif(n)
  mu    <- exp(-2.02 + 3 * x1) # 0.60 approximately
  sigma <- exp(2.23 - 2 * x2)  # 3.42 approximately
  nu    <- 2
  y <- rMOEIW(n=n, mu=mu, sigma=sigma, nu=nu)
  data.frame(y=y, x1=x1, x2=x2)
}

set.seed(123)
dat <- gendat(n=100)

mod <- gamlss(y~x1, sigma.fo=~x2, nu.fo=~1, 
              family=MOEIW, data=dat,
              control=gamlss.control(n.cyc=5000, trace=FALSE))

coef(mod, what="mu")
coef(mod, what="sigma")
exp(coef(mod, what="nu"))

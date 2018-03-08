require(RelDists)
require(gamlss)


## Estimating parameters
mu <- 0.5     # Fixing mu
sigma <- 0.1  # Fixing sigma
y <- rFWE(n=100, mu=mu, sigma=sigma)
require(gamlss)
mod <- gamlss(y ~ 1, sigma.fo= ~ 1, family=FWE)
summary(mod)
# to obtain the estimated parameters we must
# use the inverse link function
exp(coef(mod, "mu"))
exp(coef(mod, "sigma"))

## Estimating parameters using covariates
n <- 500
x1 <- runif(n)
x2 <- rpois(n, lambda=3)
mu    <- exp(-1 + 3.5 * x1)
sigma <- exp(-4.9 + 2 * x2)
y <- rFWE(n=n, mu=mu, sigma=sigma)

mod <- NULL
mod <- gamlss(y ~ x1, sigma.fo= ~ x2, family=FWE,
              control=gamlss.control(n.cyc=20, trace=FALSE),
              i.control=glim.control(cyc=50, glm.trace=FALSE))
summary(mod)

mod <- gamlss(y ~ x1, sigma.fo= ~ x2, family=FWE)

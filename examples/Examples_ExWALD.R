# Example 1
# Generating random values with
# known mu, sigma and nu

mu <- 0.20
sigma <- 70
nu <- 115

set.seed(123)
y <- rExWALD(n=100, mu, sigma, nu)

library(gamlss)
mod1 <- gamlss(y~1, family=ExWALD,
               control=gamlss.control(n.cyc=1000, trace=TRUE))

exp(coef(mod1, what="mu"))
exp(coef(mod1, what="sigma"))
exp(coef(mod1, what="nu"))

# Example 2
# Generating random values under some model

\donttest{
# A function to simulate a data set with Y ~ ExWALD
gendat <- function(n) {
  x1 <- runif(n)
  x2 <- runif(n)
  mu    <- exp(-1 + 2.8 * x1) # 1.5 approximately
  sigma <- exp( 1 - 1.2 * x2) # 1.5 approximately
  nu    <- 2
  y <- rExWALD(n=n, mu=mu, sigma=sigma, nu=nu)
  data.frame(y=y, x1=x1, x2=x2)
}

set.seed(123)
dat <- gendat(n=100)

# Fitting the model
mod2 <- gamlss(y ~ x1,
               sigma.fo = ~ x2,
               nu.fo = ~ 1,
               family = ExWALD,
               data = dat,
               control = gamlss.control(n.cyc=1000, 
                                        trace=TRUE))
summary(mod2)
}

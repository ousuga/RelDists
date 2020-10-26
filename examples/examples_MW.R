# Example 1
# Generating some random values with
# known mu, sigma and nu
y <- rMW(n=100, mu = 2, sigma = 1.5, nu = 0.2)

# Fitting the model
require(gamlss)

mod <- gamlss(y~1, sigma.fo=~1, nu.fo=~1, family= 'MW',
               control=gamlss.control(n.cyc=5000, trace=FALSE))

# Extracting the fitted values for mu, sigma and nu
# using the inverse link function
exp(coef(mod, what='mu'))
exp(coef(mod, what='sigma'))
exp(coef(mod, what='nu'))

# Example 2
# Generating random values under some model
n     <- 200
x1    <- rpois(n, lambda=2)
x2    <- runif(n)
mu    <- exp(3 -1 * x1)
sigma <- exp(2 - 2 * x2)
nu    <- 0.2
x     <- rMW(n=n, mu, sigma, nu)

mod <- gamlss(x~x1, mu.fo=~x1, sigma.fo=~x2, nu.fo=~1, family=MW,
              control=gamlss.control(n.cyc=5000, trace=FALSE))

coef(mod, what="mu")
coef(mod, what="sigma")
coef(mod, what='nu')

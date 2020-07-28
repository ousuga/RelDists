# Example 1
# Generating some random values with
# known mu and sigma
y <- rOW(n=100, mu=0.1, sigma=7, nu = 1.1)

# Fitting the model
require(gamlss)
mod <- gamlss(y~1, sigma.fo=~1, nu.fo=~1, family='OW',
              control=gamlss.control(n.cyc=5000, trace=FALSE))

# Extracting the fitted values for mu, sigma and nu
# using the inverse link function
exp(coef(mod, what='mu'))
coef(mod, what='sigma')
coef(mod, what='nu')

# Example 2
# Generating random values under some model
n <- 500
x1 <- runif(n)
x2 <- runif(n)
x3 <- runif(n)
mu <- exp(1.2 + 2 * x1)
sigma <- 2.12 + 3 * x2
nu <- exp(0.2 - x3)
x <- rOW(n=n, mu, sigma, nu)

mod <- gamlss(x~x1, sigma.fo=~x2, nu.fo=~x3, family=OW(sigma.link='identity'), 
              control=gamlss.control(n.cyc=5000, trace=FALSE))

coef(mod, what='mu')
coef(mod, what='sigma')
coef(mod, what='nu')

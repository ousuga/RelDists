# Example 1
# Generating some random values with
# known mu and sigma 
y <- rLW(n=100, mu=0, sigma=1.5)

# Fitting the model
require(gamlss)

mod <- gamlss(y~1, sigma.fo=~1, family= 'LW',
              control=gamlss.control(n.cyc=5000, trace=FALSE))

# Extracting the fitted values for mu and sigma
# using the inverse link function
coef(mod, 'mu')
exp(coef(mod, 'sigma'))

# Example 2
# Generating random values under some model
n <- 200
x1 <- runif(n, min=0.4, max=0.6)
x2 <- runif(n, min=0.4, max=0.6)
mu <- 1.5 - 3 * x1
sigma <- exp(1.4 - 2 * x2)
x <- rLW(n=n, mu, sigma)

mod <- gamlss(x~x1, sigma.fo=~x2, family=LW,
              control=gamlss.control(n.cyc=5000, trace=FALSE))

coef(mod, what="mu")
coef(mod, what="sigma")

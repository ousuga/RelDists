# Example 1
# Generating some random values with
# known mu, sigma and nu 
y <- rCS2e(n=100, mu = 0.1, sigma =1, nu=0.5)

# Fitting the model
require(gamlss)

mod <- gamlss(y~1, sigma.fo=~1, nu.fo=~1, family='CS2e',
              control=gamlss.control(n.cyc=5000, trace=FALSE))

# Extracting the fitted values for mu, sigma and nu
# using the inverse link function
exp(coef(mod, what='mu'))
exp(coef(mod, what='sigma'))
exp(coef(mod, what='nu'))

# Example 2
# Generating random values under some model
n <- 200
x1 <- runif(n, min=0.45, max=0.55)
x2 <- runif(n, min=0.4, max=0.6)
mu <- exp(0.2 - x1)
sigma <- exp(0.8 - x2)
nu <- 0.5
x <- rCS2e(n=n, mu, sigma, nu)

mod <- gamlss(x~x1, sigma.fo=~x2, nu.fo=~1,family=CS2e,
              control=gamlss.control(n.cyc=50000, trace=FALSE))

coef(mod, what="mu")
coef(mod, what="sigma")
exp(coef(mod, what="nu"))

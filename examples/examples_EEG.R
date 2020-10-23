# Generating some random values with
# known mu, sigma, nu and tau
y <- rEEG(n=100, mu = 1, sigma =1.5)

# Fitting the model
require(gamlss)

mod <- gamlss(y~1, sigma.fo=~1, family=EEG,
              control=gamlss.control(n.cyc=5000, trace=FALSE))

# Extracting the fitted values for mu, sigma, nu and tau
# using the inverse link function
exp(coef(mod, what='mu'))
exp(coef(mod, what='sigma'))

# Example 2
# Generating random values under some model
n <- 200
x1 <- runif(n, min=0.1, max=0.2)
x2 <- runif(n, min=0.1, max=0.15)
mu <- exp(0.75 - x1)
sigma <- exp(0.5 - x2)
x <- rEEG(n=n, mu, sigma)

mod <- gamlss(x~x1, sigma.fo=~x2, family=EEG,
              control=gamlss.control(n.cyc=5000, trace=FALSE))

coef(mod, what="mu")
coef(mod, what="sigma")

# Example 1
# Generating some random values with
# known mu, sigma, nu and tau
# Will not be run this example because high number is cycles
# is needed in order to get good estimates
\dontrun{
y <- rAddW(n=100, mu=1.5, sigma=0.2, nu=3, tau=0.8)

# Fitting the model
require(gamlss)

mod <- gamlss(y~1, sigma.fo=~1, nu.fo=~1, tau.fo=~1, family='AddW',
              control=gamlss.control(n.cyc=5000, trace=FALSE))

# Extracting the fitted values for mu, sigma, nu and tau
# using the inverse link function
exp(coef(mod, what='mu'))
exp(coef(mod, what='sigma'))
exp(coef(mod, what='nu'))
exp(coef(mod, what='tau'))
}

# Example 2
# Generating random values under some model
# Will not be run this example because high number is cycles
# is needed in order to get good estimates
\dontrun{
n <- 200
x1 <- runif(n, min=0.4, max=0.6)
x2 <- runif(n, min=0.4, max=0.6)
mu <- exp(1.67 + -3 * x1)
sigma <- exp(0.69 - 2 * x2)
nu <- 3
tau <- 0.8
x <- rAddW(n=n, mu, sigma, nu, tau)

mod <- gamlss(x~x1, sigma.fo=~x2, nu.fo=~1, tau.fo=~1, family=AddW,
              control=gamlss.control(n.cyc=5000, trace=FALSE))

coef(mod, what="mu")
coef(mod, what="sigma")
exp(coef(mod, what="nu"))
exp(coef(mod, what="tau"))
}

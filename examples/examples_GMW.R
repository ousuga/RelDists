# Example 1
# Generating some random values with
# known mu, sigma, nu and tau
y <- rGMW(n=100, mu=2, sigma=0.5, nu=2, tau=1.5)

# Fitting the model
require(gamlss)

mod <- gamlss(y~1, sigma.fo=~1, nu.fo=~1, tau.fo=~ 1, family='GMW',
              control=gamlss.control(n.cyc=5000, trace=FALSE))

# Extracting the fitted values for mu, sigma and nu
# using the inverse link function
exp(coef(mod, what='mu'))
exp(coef(mod, what='sigma'))
(coef(mod, what='nu'))^2
(coef(mod, what='tau'))^2

# Example 2
# Generating random values under some model
\dontrun{
n <- 1000
x1 <- runif(n)
x2 <- runif(n)
mu <- exp(2 + -3 * x1)
sigma <- exp(3 - 2 * x2)
nu <- 2
tau <- 1.5
x <- rGMW(n=n, mu, sigma, nu, tau)

mod <- gamlss(x~x1, sigma.fo=~x2, nu.fo=~1, tau.fo=~ 1, family="GMW", 
              control=gamlss.control(n.cyc=5000, trace=FALSE))

coef(mod, what="mu")
coef(mod, what="sigma")
coef(mod, what="nu")^2
coef(mod, what="tau")^2
}

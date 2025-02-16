# Example 1
# Generating some random values with
# known mu and sigma
y <- rIW(n=100, mu=1, sigma=2)

# Fitting the model
require(gamlss)
mod <- gamlss(y~1, mu.fo=~1, sigma.fo=~1, family="IW")

# Extracting the fitted values for mu, sigma and nu
# using the inverse link function
exp(coef(mod, what="mu"))
exp(coef(mod, what="sigma"))

# Example 2
# Generating random values under some model
n <- 100
x1 <- runif(n)
x2 <- runif(n)
mu <- exp(2 + -1 * x1)
sigma <- exp(2 - 2 * x2)
y <- rIW(n=n, mu=mu, sigma=sigma)

mod <- gamlss(y~x1, mu.fo=~1, sigma.fo=~x2, family=IW)

coef(mod, what="mu")
coef(mod, what="sigma")

# Example 3
# Using the dataset from Kundu and Howlader (2010) Bayesian inference and 
# prediction of the inverse Weibull distribution for Type-II censored data

y <- c(12, 15, 22, 24, 24, 32, 32, 33, 34, 38, 38, 43, 44, 48, 52, 
       53, 54, 54, 55, 56, 57, 58, 58, 59, 60, 60, 60, 60, 61, 62, 
       63, 65, 65, 67, 68, 70, 70, 72, 73, 75, 76, 76, 81, 83, 84, 
       85, 87, 91, 95, 96, 98, 99, 109, 110, 121, 127, 129, 131, 
       143, 146, 146, 175, 175, 211, 233, 258, 258, 263, 297, 341, 
       341, 376)

y <- y / 1000

mod <- gamlss(y~1, mu.fo=~1, sigma.fo=~1, family="IW",
              control=gamlss.control(n.cyc=5000, trace=FALSE))
exp(coef(mod, what="mu"))
exp(coef(mod, what="sigma"))



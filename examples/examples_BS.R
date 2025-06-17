# Example 1
# Generating some random values with
# known mu and sigma
y <- rBS(n=100, mu=0.75, sigma=1.3)

# Fitting the model
require(gamlss)
mod1 <- gamlss(y~1, sigma.fo=~1, family=BS)

# Extracting the fitted values for mu and sigma
# using the inverse link function
exp(coef(mod1, what="mu"))
exp(coef(mod1, what="sigma"))


# Example 2
# Generating random values for a regression model

# A function to simulate a data set with Y ~ BS
gendat <- function(n) {
  x1 <- runif(n)
  x2 <- runif(n)
  mu <- exp(1.45 - 3 * x1)
  sigma <- exp(2 - 1.5 * x2)
  y <- rBS(n=n, mu=mu, sigma=sigma)
  data.frame(y=y, x1=x1, x2=x2)
}

set.seed(123)
dat <- gendat(n=300)

mod2 <- gamlss(y~x1, sigma.fo=~x2, 
               family=BS, data=dat)

summary(mod2)

# Example 3

# Example 1
# Generating some random values with
# known mu and sigma
y <- rBS3(n=50, mu=2, sigma=0.2)

# Fitting the model
require(gamlss)
mod1 <- gamlss(y~1, sigma.fo=~1, family=BS3)

# Extracting the fitted values for mu and sigma
# using the inverse link function
exp(coef(mod1, what="mu"))
exp(coef(mod1, what="sigma"))

# Example 2
# Generating random values for a regression model

# A function to simulate a data set with Y ~ BS3
\dontrun{
gendat <- function(n) {
  x1 <- runif(n)
  x2 <- runif(n)
  mu <- exp(1.45 - 3 * x1)
  inv_logit <- function(x) 1 / (1 + exp(-x))
  sigma <- inv_logit(2 - 1.5 * x2)
  y <- rBS3(n=n, mu=mu, sigma=sigma)
  data.frame(y=y, x1=x1, x2=x2)
}

set.seed(1234)
dat <- gendat(n=100)

mod2 <- gamlss(y~x1, sigma.fo=~x2, 
               family=BS3, data=dat,
               control=gamlss.control(n.cyc=100))

summary(mod2)
}

# Example 3
# The response variable is the ratio between the average
# rent per acre planted with alfalfa and the corresponding 
# average rent for other agricultural uses. The density of
# dairy cows (X2, number per square mile) is the explanatory variable. 
library(alr4)
data("landrent")

landrent$ratio <- landrent$Y / landrent$X1

with(landrent, plot(x=X2, y=ratio))

mod3 <- gamlss(ratio~X2, sigma.fo=~X2, 
               data=landrent, family=BS3)

summary(mod3)
logLik(mod3)


# Example 1
# Generating some random values with
# known mu and sigma
y <- rEXL(n=300, mu=0.75, sigma=1.3)

# Fitting the model
require(gamlss)
mod1 <- gamlss(y~1, sigma.fo=~1, family=EXL,
              control=gamlss.control(n.cyc=5000, trace=FALSE))

# Extracting the fitted values for mu and sigma
# using the inverse link function
exp(coef(mod1, what="mu"))
exp(coef(mod1, what="sigma"))


# Example 2
# Generating random values for a regression model

# A function to simulate a data set with Y ~ EXL
gendat <- function(n) {
  x1 <- runif(n)
  x2 <- runif(n)
  mu <- exp(1.45 - 3 * x1)
  sigma <- exp(2 - 1.5 * x2)
  y <- rEXL(n=n, mu=mu, sigma=sigma)
  data.frame(y=y, x1=x1, x2)
}

set.seed(1234)
dat <- gendat(n=100)

mod2 <- gamlss(y~x1, sigma.fo=~x2, 
               family=EXL, data=dat,
               control=gamlss.control(n.cyc=5000, trace=FALSE))

summary(mod2)

# Example 3
# Mortality rate due to COVID-19 for 30 days (31st March to April 30, 2020)
# recorded for the Netherlands.
# Taken from Alomair et al. (2024) page 12.

x <- c(14.918, 10.656, 12.274, 10.289, 10.832, 7.099, 5.928, 13.211, 
       7.968, 7.584, 5.555, 6.027, 4.097, 3.611, 4.960, 7.498, 6.940, 
       5.307, 5.048, 2.857, 2.254, 5.431, 4.462, 3.883,
       3.461, 3.647, 1.974, 1.273, 1.416, 4.235)

mod3 <- gamlss(x~1, sigma.fo=~1, family=EXL,
               control=gamlss.control(n.cyc=5000, trace=FALSE))

# Extracting the fitted values for mu and sigma
# using the inverse link function
exp(coef(mod3, what="mu"))
exp(coef(mod3, what="sigma"))

# Replicating figure 4 from Alomair et al. (2024)
# Hist and estimated pdf
hist(x, freq=FALSE)
curve(dEXL(x, mu=0.4089915, sigma=2.710467), add=TRUE, col="tomato", lwd=2)
# Empirical cdf and estimated ecdf
plot(ecdf(x))
curve(pEXL(x, mu=0.4089915, sigma=2.710467), add=TRUE, col="tomato", lwd=2)
# QQplot
qqplot(x, rEXL(n=30, mu=0.4089915, sigma=2.710467), col="tomato")
qqline(x, distribution=function(p) qEXL(p, mu=0.4089915, sigma=2.710467))


# Example 4
# Precipitation in inches
# Taken from Alomair et al. (2024) page 13.

# Manuel
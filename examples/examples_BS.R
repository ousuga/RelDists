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
# Fatigue life (T) measures in cycles (×10−3) of n 101
# aluminum coupons (specimens) of type 6061-T6.
# Taken from Leiva et al. (2006) page 37.
# https://journal.r-project.org/articles/RN-2006-033/RN-2006-033.pdf

y <- c(70, 90, 96, 97, 99, 100, 103, 104,
       104, 105, 107, 108, 108, 108, 109, 109,
       112, 112, 113, 114, 114, 114, 116, 119,
       120, 120, 120, 121, 121, 123, 124, 124,
       124, 124, 124, 128, 128, 129, 129, 130,
       130, 130, 131, 131, 131, 131, 131, 132,
       132, 132, 133, 134, 134, 134, 134, 134,
       136, 136, 137, 138, 138, 138, 139, 139,
       141, 141, 142, 142, 142, 142, 142, 142,
       144, 144, 145, 146, 148, 148, 149, 151,
       151, 152, 155, 156, 157, 157, 157, 157,
       158, 159, 162, 163, 163, 164, 166, 166,
       168, 170, 174, 196, 212)

mod3 <- gamlss(y~1, sigma.fo=~1, family=BS)

# Extracting the fitted values for mu and sigma
# using the inverse link function
exp(coef(mod3, what="mu"))
exp(coef(mod3, what="sigma"))

# Example 4
# Aggregate payments by the insurer
# in thousand Skr (Swedish currency).
# Taken from Balakrishnan and Kundu (2019) page 65.
# https://onlinelibrary.wiley.com/doi/abs/10.1002/asmb.2348

y <- c(5014, 5855, 6486, 6540, 6656, 6656, 7212, 7541, 7558, 
       7797, 8546, 9345, 11762, 12478, 13624, 14451,
       14940, 14963, 15092, 16203, 16229, 16730, 18027, 
       18343, 19365, 21782, 24248, 29069, 34267, 38993)

y <- y/10000

mod4 <- gamlss(y~1, sigma.fo=~1, family=BS)

# Extracting the fitted values for mu and sigma
# using the inverse link function
exp(coef(mod4, what="mu"))
exp(coef(mod4, what="sigma"))




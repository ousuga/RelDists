# Example 1
# Generating some random values with
# known mu, sigma and nu
y <- rOW(n=200, mu=0.2, sigma=4, nu=0.05)

# Custom search region
myvalues <- list(sigma="all(sigma > 1)",
                 nu="all(nu < 1) & all(nu < 1)")

require(gamlss)
mod <- gamlss(y~1, sigma.fo=~1, nu.fo=~1, 
              sigma.start=2, nu.start=0.1,
              control=gamlss.control(n.cyc=300, trace=FALSE),
              family=myOW_region(valid.values=myvalues))

exp(coef(mod, what='mu'))
exp(coef(mod, what='sigma'))
exp(coef(mod, what='nu'))

# Example 2
# Same example using another link function
mod <- gamlss(y~1, sigma.fo=~1, nu.fo=~1, 
              sigma.start=2, nu.start=0.1,
              control=gamlss.control(n.cyc=300, trace=FALSE),
              family=myOW_region(family=OW(sigma.link='identity'),
                                 valid.values=myvalues))

exp(coef(mod, what='mu'))
coef(mod, what='sigma')
exp(coef(mod, what='nu'))

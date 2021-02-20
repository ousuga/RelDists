# Example 1
# Generating some random values with
# known mu, sigma and nu
y <- rOW(n=200, mu=0.2, sigma=4, nu=0.05)

# Custom search region
myvalues <- list(sigma="all(sigma > 1)",
                 nu="all(nu < 1) & all(nu < 1)")

my_initial_guess <- initValuesOW_TTT(formula=y~1)
summary(my_initial_guess)

# OW family modified with 'myOW_region'
require(gamlss)
myOW <- myOW_region(valid.values=myvalues, initVal=my_initial_guess)
mod1 <- gamlss(y~1, sigma.fo=~1, nu.fo=~1, 
               sigma.start=2, nu.start=0.1,
               control=gamlss.control(n.cyc=300, trace=FALSE),
               family=myOW)

exp(coef(mod1, what='mu'))
exp(coef(mod1, what='sigma'))
exp(coef(mod1, what='nu'))

# Example 2
# Same example using another link function and using 'myOW_region'
# in the argument 'family'
mod2 <- gamlss(y~1, sigma.fo=~1, nu.fo=~1, 
               sigma.start=2, nu.start=0.1,
               control=gamlss.control(n.cyc=300, trace=FALSE),
               family=myOW_region(family=OW(sigma.link='identity'),
                                  valid.values=myvalues,
                                  initVal=my_initial_guess))

exp(coef(mod2, what='mu'))
coef(mod2, what='sigma')
exp(coef(mod2, what='nu'))

# Example 3
# Fitting an OW model with covariates
set.seed(321)
x <- c(rep(0,30), rep(1,30))
nu <- 2
sigma <- 0.7 + 2*x
probs <- runif(n=length(x))
y2 <- qOW(p=probs, mu=0.005, sigma=sigma, nu=nu)
x <- as.factor(x)

myoptions <- list(loess.options(), loess.options(span = 0.75))
my_initial_val <- initValuesOW_TTT(y2 ~ x, local_reg = myoptions)
summary(my_initial_val)
plot(my_initial_val, legend_options=list(pos=1.03),
     par_plot=list(mar=c(3.7,3.7,1,2.8), mgp=c(2.5,1,0)))

mod3 <- gamlss(y2~1, sigma.fo=~x, nu.fo=~1, 
               sigma.start=2, nu.start=0.1,
               control=gamlss.control(n.cyc=300, trace=FALSE),
               family=myOW_region(family=OW,
                                  initVal=my_initial_val, level=2))
exp(coef(mod3, what='mu'))
exp(coef(mod3, what='sigma'))
exp(coef(mod3, what='nu'))

# Random data generation (OW distributed)
y <- rOW(n=1000, mu=0.05, sigma=0.6, nu=2)

# Initial values with TTT plot
iv <- initValuesOW_TTT(formula = y ~ 1)
summary(iv)

# This data is from unimodal hazard
# See TTT estimate from sample
plot(iv, legend_options=list(pos=1.03))

# See the true hazard
curve(hOW(x, mu=0.05, sigma=0.6, nu=2), to=100, lwd=3, ylab="h(x)")

# Finally, we fit the model
require(gamlss)
con.out <-gamlss.control(n.cyc = 300, trace = FALSE)
con.in <- glim.control(cyc = 500)

sigma.start <- param.startOW("sigma", iv)
nu.start <- param.startOW("nu", iv)

mod <- gamlss(y~1, sigma.fo=~1, nu.fo=~1, control=con.out, i.control=con.in,
              family=myOW_region(OW(sigma.link="identity", nu.link="identity"),
                                 valid.values="auto", iv),
              sigma.start=sigma.start, nu.start=nu.start)

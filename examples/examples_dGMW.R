## The probability density function
curve(dGMW(x, mu=2, sigma=0.5, nu=2, tau=1.5), from=0, to=0.8,
      ylim=c(0, 2), col="red", las=1, ylab="f(x)") 
 
## The cumulative distribution and the Reliability function
par(mfrow=c(1, 2))
curve(pGMW(x, mu=2, sigma=0.5, nu=2, tau=1.5),
      from=0, to=1.2, col="red", las=1, ylab="F(x)")
curve(pGMW(x, mu=2, sigma=0.5, nu=2, tau=1.5, lower.tail=FALSE),
      from=0, to=1.2, col="red", las=1, ylab="R(x)")

## The quantile function
p <- seq(from=0, to=0.99999, length.out=100)
plot(x=qGMW(p, mu=2, sigma=0.5, nu=2, tau=1.5), y=p, xlab="Quantile",
    las=1, ylab="Probability")
curve(pGMW(x, mu=2, sigma=0.5, nu=2, tau=1.5), from=0, add=TRUE, col="red")
 
## The random function
hist(rGMW(n=1000, mu=2, sigma=0.5, nu=2,tau=1.5), freq=FALSE,
    xlab="x", main ="", las=1)
curve(dGMW(x, mu=2, sigma=0.5, nu=2, tau=1.5),  from=0, add=TRUE, col="red")
 
## The Hazard function
par(mfrow=c(1,1))
curve(hGMW(x, mu=2, sigma=0.5, nu=2, tau=1.5), from=0, to=1, ylim=c(0, 16),
      col="red", ylab="The Hazard function", las=1)
      
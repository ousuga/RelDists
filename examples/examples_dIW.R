## The probability density function
curve(dIW(x, mu=5, sigma=2.5), from=0, to=10,
      ylim=c(0, 0.55), col="red", las=1, ylab="f(x)")
#'
## The cumulative distribution and the Reliability function
par(mfrow=c(1, 2))
curve(pIW(x, mu=5, sigma=2.5),
      from=0, to=10,  col="red", las=1, ylab="F(x)")
curve(pIW(x, mu=5, sigma=2.5, lower.tail=FALSE),
      from=0, to=10, col="red", las=1, ylab="R(x)")
            
## The quantile function
p <- seq(from=0, to=0.99999, length.out=100)
plot(x=qIW(p, mu=5, sigma=2.5), y=p, xlab="Quantile",
  las=1, ylab="Probability")
curve(pIW(x, mu=5, sigma=2.5), from=0, add=TRUE, col="red")
  
## The random function
hist(rIW(n=10000, mu=5, sigma=2.5), freq=FALSE, xlim=c(0,60),
  xlab="x", las=1, main="")
curve(dIW(x, mu=5, sigma=2.5), from=0, add=TRUE, col="red")

## The Hazard function
par(mfrow=c(1,1))
curve(hIW(x, mu=5, sigma=2.5), from=0, to=15, ylim=c(0, 0.9),
   col="red", ylab="Hazard function", las=1)

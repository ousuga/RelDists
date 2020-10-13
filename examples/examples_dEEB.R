## The probability density function
par(mfrow = c(1, 2))
curve(dEEB(x, mu=1, sigma=1, nu=0.5, size=10), from=0, to=1.5, ylim=c(0, 1.5), 
      col="red", las=1, ylab="f(x)")
curve(dEEB(x, mu=1, sigma=1, nu=5, size=10), from=0, to=1.5, ylim=c(0, 3), 
      col="red", las=1, ylab="f(x)")
      
## The cumulative distribution and the Reliability function
par(mfrow = c(1, 2))
curve(pEEB(x, mu=1, sigma=1, nu=5, size=10), from=0, to=4, ylim=c(0, 1), 
      col="red", las = 1, ylab = "F(x)")
curve(pEEB(x, mu=1, sigma=1, nu=5, size=10, lower.tail=FALSE), from=0, to=4,  
      ylim=c(0, 1), col="red", las=1, ylab="R(x)")

## The quantile function
p <- seq(from=0, to=0.998, length.out=100)
plot(x=qEEB(p, mu=1, sigma=1, nu=5, size=10), y=p, xlab="Quantile", 
    las=1, ylab="Probability")
curve(pEEB(x, mu=1, sigma=1, nu=5, size=10), from=0, add=TRUE, col="red")

## The random function
hist(rEEB(n=10000, mu=1, sigma=1, nu=5, size=10), freq=FALSE, 
     ylim = c(0,1.5),xlab="x", las=1, main="")
curve(dEEB(x, mu=1, sigma=1, nu=5, size=10), from=0, ylim=c(0, 3), 
      add=TRUE, col="red")

## The Hazard function
par(mfrow=c(1,2))
curve(hEEB(x, mu=1, sigma=1, nu=0.5, size=10), from=0, to=4, ylim=c(0, 30), 
      col="red", ylab="Hazard function", las=1)
curve(hEEB(x, mu=1, sigma=1, nu=5, size=10), from=0, to=6, ylim=c(0, 10), 
      col="red", ylab="Hazard function", las=1)


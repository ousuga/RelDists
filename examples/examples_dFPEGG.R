## The probability density function
curve(dFPEGG(x, mu=0.1, sigma=0.8, nu=10, tau=1.5), from=0.000001, to=1.5, ylim=c(0, 2.5),
      col="red", las=1, ylab="f(x)")

## The cumulative distribution and the Reliability function
par(mfrow=c(1, 2))
curve(pFPEGG(x, mu=0.1, sigma=0.8, nu=10, tau=1.5),
      from=0.000001, to=1.5, col="red", las=1, ylab="F(x)")
curve(pFPEGG(x, mu=0.1, sigma=0.8, nu=10, tau=1.5, lower.tail=FALSE),
      from=0.000001, to=1.5, col="red", las=1, ylab="R(x)")

## The quantile function
p <- seq(from=0, to=0.99999, length.out=100)
plot(x=qFPEGG(p, mu=0.1, sigma=0.8, nu=10, tau=1.5), y=p, xlab="Quantile",
     las=1, ylab="Probability")
curve(pFPEGG(x, mu=0.1, sigma=0.8, nu=10, tau=1.5), 
      from=0.00001, add=TRUE, col="red")

## The random function
hist(rFPEGG(n=100, mu=0.1, sigma=0.8, nu=10, tau=1.5), freq=FALSE,
     xlab="x", las=1, main="")
curve(dFPEGG(x, mu=0.1, sigma=0.8, nu=10, tau=1.5),
      from=0.0001, to=2, add=TRUE, col="red")

## The Hazard function
curve(hFPEGG(x,  mu=0.1, sigma=0.8, nu=10, tau=1.5), from=0.0001, to=1.5,
      col="red", ylab="Hazard function", las=1)


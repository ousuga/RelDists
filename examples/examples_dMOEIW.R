## The probability density function
curve(dMOEIW(x, mu=0.6, sigma=1.7, nu=0.3), from=0, to=2,
      col="red", ylab="f(x)", las=1)

## The cumulative distribution and the Reliability function
par(mfrow=c(1, 2))
curve(pMOEIW(x, mu=0.6, sigma=1.7, nu=0.3),
      from=0.0001, to=2, col="red", las=1, ylab="F(x)")
curve(pMOEIW(x, mu=0.6, sigma=1.7, nu=0.3, lower.tail=FALSE),
      from=0.0001, to=2, col="red", las=1, ylab="R(x)")

## The quantile function
p <- seq(from=0, to=0.99999, length.out=100)
plot(x=qMOEIW(p, mu=0.6, sigma=1.7, nu=0.3), y=p, xlab="Quantile",
     las=1, ylab="Probability")
curve(pMOEIW(x, mu=0.6, sigma=1.7, nu=0.3),
      from=0, add=TRUE, col="red")

## The random function
hist(rMOEIW(n=1000, mu=0.6, sigma=1.7, nu=0.3), freq=FALSE,
     xlab="x", las=1, main="")
curve(dMOEIW(x, mu=0.6, sigma=1.7, nu=0.3),
      from=0.001, to=4, add=TRUE, col="red")

## The Hazard function
curve(hMOEIW(x, mu=0.5, sigma=0.7, nu=1), from=0.001, to=3,
      col="red", ylab="Hazard function", las=1)

## The probability density function
curve(dGIW(x, mu=3, sigma=5, nu=0.5), from=0.001, to=8,
      col="red", ylab="f(x)", las=1)

## The cumulative distribution and the Reliability function
par(mfrow=c(1, 2))
curve(pGIW(x, mu=3, sigma=5, nu=0.5),
      from=0.0001, to=14, col="red", las=1, ylab="F(x)")
curve(pGIW(x, mu=3, sigma=5, nu=0.5, lower.tail=FALSE),
      from=0.0001, to=14, col="red", las=1, ylab="R(x)")

## The quantile function
p <- seq(from=0, to=0.99999, length.out=100)
plot(x=qGIW(p, mu=3, sigma=5, nu=0.5), y=p, xlab="Quantile",
     las=1, ylab="Probability")
curve(pGIW(x, mu=3, sigma=5, nu=0.5),
      from=0, add=TRUE, col="red")

## The random function
hist(rGIW(n=1000, mu=3, sigma=5, nu=0.5), freq=FALSE,
     xlab="x", ylim=c(0, 0.8), las=1, main="")
curve(dGIW(x, mu=3, sigma=5, nu=0.5),
      from=0.001, to=14, add=TRUE, col="red")

## The Hazard function
curve(hGIW(x, mu=3, sigma=5, nu=0.5), from=0.001, to=30,
      col="red", ylab="Hazard function", las=1)
#


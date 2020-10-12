## The probability density function
curve(dWP(x, mu=1.5, sigma=0.5, nu=10), from=0.0001, to=2,
      col="red", las=1, ylab="f(x)")

## The cumulative distribution and the Reliability function
par(mfrow=c(1, 2))
curve(pWP(x, mu=1.5, sigma=0.5, nu=10),
      from=0.0001, to=2, col="red", las=1, ylab="F(x)")
curve(pWP(x, mu=1.5, sigma=0.5, nu=10, lower.tail=FALSE),
      from=0.0001, to=2, col="red", las=1, ylab="R(x)")

## The quantile function
p <- seq(from=0, to=0.99999, length.out=100)
plot(x=qWP(p, mu=1.5, sigma=0.5, nu=10), y=p, xlab="Quantile",
     las=1, ylab="Probability")
curve(pWP(x, mu=1.5, sigma=0.5, nu=10),
      from=0, add=TRUE, col="red")

## The random function
hist(rWP(n=10000, mu=1.5, sigma=0.5, nu=10), freq=FALSE,
     xlab="x", ylim=c(0, 2.2), las=1, main="")
curve(dWP(x, mu=1.5, sigma=0.5, nu=10),
      from=0.001, to=4, add=TRUE, col="red")

## The Hazard function
curve(hWP(x, mu=1.5, sigma=0.5, nu=10), from=0.001, to=5,
      col="red", ylab="Hazard function", las=1)

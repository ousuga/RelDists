## The probability density function
curve(dQXGP(x, mu=0.5, sigma=1, nu=1), from=0.1, to=8,
      ylim=c(0, 0.6), col="red", las=1, ylab="f(x)")

## The cumulative distribution and the Reliability function
par(mfrow=c(1, 2))
curve(pQXGP(x, mu=0.5, sigma=1, nu=1),
      from=0.1, to=8,  col="red", las=1, ylab="F(x)")
curve(pQXGP(x,  mu=0.5, sigma=1, nu=1, lower.tail=FALSE),
      from=0.1, to=8, col="red", las=1, ylab="R(x)")

## The quantile function
p <- seq(from=0, to=0.99999, length.out=100)
plot(x=qQXGP(p, mu=0.5, sigma=1, nu=1), y=p, xlab="Quantile",
     las=1, ylab="Probability")
curve(pQXGP(x, mu=0.5, sigma=1, nu=1),
      from=0.1, add=TRUE, col="red")
      
## The random function
hist(rQXGP(n=1000, mu=0.5, sigma=1, nu=1), freq=FALSE,
     xlab="x", ylim=c(0, 0.4), las=1, main="", xlim=c(0, 15))
curve(dQXGP(x, mu=0.5, sigma=1, nu=1),
      from=0.001, to=500, add=TRUE, col="red")

## The Hazard function
curve(hQXGP(x, mu=0.5, sigma=1, nu=1), from=0.01, to=3,
      col="red", ylab="Hazard function", las=1)

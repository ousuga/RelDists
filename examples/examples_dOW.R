## The probability density function
curve(dOW(x, mu=2, sigma=3, nu=0.2), from=0, to=4, ylim=c(0, 2), 
      col="red", las=1, ylab="f(x)")

## The cumulative distribution and the Reliability function
par(mfrow = c(1, 2))
curve(pOW(x, mu=2, sigma=3, nu=0.2), from=0, to=4, ylim=c(0, 1), 
      col="red", las=1, ylab="f(x)")
curve(pOW(x, mu=2, sigma=3, nu=0.2, lower.tail=FALSE), from=0, 
      to=4,  ylim=c(0, 1), col="red", las=1, 
      ylab = "S(x)")

## The quantile function
p <- seq(from=0, to=0.998, length.out=100)
plot(x = qOW(p, mu=2, sigma=3, nu=0.2), y=p, xlab="Quantile", las=1, 
     ylab="Probability")
curve(pOW(x, mu=2, sigma=3, nu=0.2), from=0, add=TRUE, col="red")

## The random function
hist(rOW(n=10000, mu=2, sigma = 3, nu = 0.2), freq=FALSE, ylim = c(0, 2),
     xlab = "x", las = 1, main = "")
curve(dOW(x, mu=2, sigma=3, nu=0.2),  from=0, ylim=c(0, 2), add=TRUE,
      col = "red")

## The Hazard function
par(mfrow=c(1,1))
curve(hOW(x, mu = 2, sigma=3, nu=0.2), from=0, to=2.5, ylim=c(0, 3),
      col="red", ylab="Hazard function", las=1)

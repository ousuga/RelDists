old_par <- par(mfrow = c(1, 1)) # save previous graphical parameters

## The probability density function 
curve(dMW(x, mu=2, sigma=1.5, nu=0.2), from=0, to=2,
  ylim=c(0, 1.5), col="red", las=1, ylab="f(x)")

## The cumulative distribution and the Reliability function
par(mfrow = c(1, 2))
curve(pMW(x, mu=2, sigma=1.5, nu=0.2), from=0, to=2,
 col = "red", las=1, ylab="F(x)")
curve(pMW(x, mu=2, sigma=1.5, nu=0.2, lower.tail = FALSE), 
from=0, to=2, col="red", las=1, ylab ="R(x)")

## The quantile function
p <- seq(from=0, to=0.9999, length.out=100)
plot(x=qMW(p, mu=2, sigma=1.5, nu=0.2), y=p, xlab="Quantile",
 las=1, ylab="Probability")
curve(pMW(x, mu=2, sigma=1.5, nu=0.2), from=0, add=TRUE, col="red")

## The random function
hist(rMW(n=1000, mu=2, sigma=1.5, nu=0.2), freq=FALSE,
 xlab="x", las=1, main="")
curve(dMW(x, mu=2, sigma=1.5, nu=0.2), from=0, add=TRUE, col="red")

## The Hazard function
par(mfrow=c(1,1))
curve(hMW(x, mu=2, sigma=1.5, nu=0.2), from=0, to=1.5, ylim=c(0, 5),
 col="red", las=1, ylab="H(x)", las=1)

par(old_par) # restore previous graphical parameters
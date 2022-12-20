old_par <- par(mfrow = c(1, 1)) # save previous graphical parameters

## The probability density function
curve(dPL(x, mu=1.5, sigma=0.2), from=0.1, to=10,
      col="red", las=1, ylab="f(x)")

## The cumulative distribution and the Reliability function
par(mfrow=c(1, 2))
curve(pPL(x, mu=1.5, sigma=0.2),
      from=0.1, to=10, col="red", las=1, ylab="F(x)")
curve(pPL(x, mu=1.5, sigma=0.2, lower.tail=FALSE),
      from=0.1, to=10, col="red", las=1, ylab="R(x)")

## The quantile function
p <- seq(from=0, to=0.99999, length.out=100)
plot(x=qPL(p, mu=1.5, sigma=0.2), y=p, xlab="Quantile",
     las=1, ylab="Probability")
curve(pPL(x, mu=1.5, sigma=0.2), from=0.1, add=TRUE, col="red")

## The random function
hist(rPL(n=1000, mu=1.5, sigma=0.2), freq=FALSE,
     xlab="x", las=1, main="")
curve(dPL(x, mu=1.5, sigma=0.2), from=0.1, to=15, add=TRUE, col="red")

## The Hazard function
par(mfrow=c(1,1))
curve(hPL(x, mu=1.5, sigma=0.2), from=0.1, to=15,
      col="red", ylab="Hazard function", las=1)

par(old_par) # restore previous graphical parameters
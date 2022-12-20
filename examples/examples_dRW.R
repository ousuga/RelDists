old_par <- par(mfrow = c(1, 1)) # save previous graphical parameters

## The probability density function
curve(dRW(x, mu=1, sigma=1), from=-5, to=-0.01,
      col="red", las=1, ylab="f(x)")

## The cumulative distribution and the Reliability function
par(mfrow=c(1, 2))
curve(pRW(x, mu=1, sigma=1),
      from=-5, to=-0.01, col="red", las=1, ylab="F(x)")
curve(pRW(x, mu=1, sigma=1, lower.tail=FALSE),
      from=-5, to=-0.01, col="red", las=1, ylab="R(x)")

## The quantile function
p <- seq(from=0, to=0.99999, length.out=100)
plot(x=qRW(p, mu=1, sigma=1), y=p, xlab="Quantile",
     las=1, ylab="Probability")
curve(pRW(x, mu=1, sigma=1), from=-5, add=TRUE, col="red")

## The random function
hist(rRW(n=10000, mu=1, sigma=1), freq=FALSE,
     xlab="x", las=1, main="")
curve(dRW(x, mu=1, sigma=1), from=-5, to=-0.01, add=TRUE, col="red")

## The Hazard function
par(mfrow=c(1,1))
curve(hRW(x, mu=1, sigma=1), from=-0.3, to=-0.01,
      col="red", ylab="Hazard function", las=1)

par(old_par) # restore previous graphical parameters
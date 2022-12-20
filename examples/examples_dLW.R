old_par <- par(mfrow = c(1, 1)) # save previous graphical parameters

## The probability density function
curve(dLW(x, mu=0, sigma=1.5), from=-8, to=5,
      col="red", las=1, ylab="f(x)")

## The cumulative distribution and the Reliability function
par(mfrow=c(1, 2))
curve(pLW(x, mu=0, sigma=1.5),
      from=-8, to=5,  col="red", las=1, ylab="F(x)")
curve(pLW(x, mu=0, sigma=1.5, lower.tail=FALSE),
      from=-8, to=5, col="red", las=1, ylab="R(x)")

## The quantile function
p <- seq(from=0, to=0.99999, length.out=100)
plot(x=qLW(p, mu=0, sigma=1.5), y=p, xlab="Quantile",
     las=1, ylab="Probability")
curve(pLW(x, mu=0, sigma=1.5), from=-8, to=5, add=TRUE, col="red")

## The random function
hist(rLW(n=10000, mu=0, sigma=1.5), freq=FALSE,
     xlab="x", las=1, main="")
curve(dLW(x, mu=0, sigma=1.5), from=-15, to=6, add=TRUE, col="red")

## The Hazard function
par(mfrow=c(1,1))
curve(hLW(x, mu=0, sigma=1.5), from=-8, to=7,
      col="red", ylab="Hazard function", las=1)

par(old_par) # restore previous graphical parameters
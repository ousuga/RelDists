old_par <- par(mfrow = c(1, 1)) # save previous graphical parameters

## The probability density function
curve(dExW(x, mu=0.3, sigma=2, nu=0.05), from=0.0001, to=2,
      col="red", las=1, ylab="f(x)")

## The cumulative distribution and the Reliability function
par(mfrow=c(1, 2))
curve(pExW(x, mu=0.3, sigma=2, nu=0.05),
      from=0.0001, to=2, col="red", las=1, ylab="F(x)")
curve(pExW(x, mu=0.3, sigma=2, nu=0.05, lower.tail=FALSE),
      from=0.0001, to=2, col="red", las=1, ylab="R(x)")

## The quantile function
p <- seq(from=0, to=0.99999, length.out=100)
plot(x=qExW(p, mu=0.3, sigma=2, nu=0.05), y=p, xlab="Quantile",
     las=1, ylab="Probability")
curve(pExW(x, mu=0.3, sigma=2, nu=0.05), 
      from=0, add=TRUE, col="red")

## The random function
hist(rExW(n=10000, mu=0.3, sigma=2, nu=0.05), freq=FALSE,
     xlab="x", ylim=c(0, 2), las=1, main="")
curve(dExW(x, mu=0.3, sigma=2, nu=0.05),
      from=0.001, to=4, add=TRUE, col="red")

## The Hazard function
curve(hExW(x, mu=0.3, sigma=2, nu=0.05), from=0.001, to=4,
      col="red", ylab="Hazard function", las=1)

par(old_par) # restore previous graphical parameters
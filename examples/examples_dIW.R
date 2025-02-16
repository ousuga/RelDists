# The probability density function
curve(dIW(x, mu=1, sigma=2), from=0, to=10,
      col="red", las=1, ylab="f(x)")

# The cumulative distribution and the Reliability function
par(mfrow=c(1, 2))
curve(pIW(x, mu=1, sigma=2),
      from=0, to=10,  col="red", las=1, ylab="F(x)")
curve(pIW(x, mu=1, sigma=2, lower.tail=FALSE),
      from=0, to=10, col="red", las=1, ylab="R(x)")
            
# The quantile function
p <- seq(from=0, to=0.99, length.out=100)
plot(x=qIW(p, mu=1, sigma=2), y=p, xlab="Quantile",
  las=1, ylab="Probability")
curve(pIW(x, mu=1, sigma=2), from=0, add=TRUE, col="red")
  
# The random function
hist(rIW(n=1000, mu=1, sigma=2), freq=FALSE, xlim=c(0, 40),
  xlab="x", las=1, main="")
curve(dIW(x, mu=1, sigma=2), from=0, add=TRUE, col="red")

# The Hazard function
par(mfrow=c(1, 1))
curve(hIW(x, mu=1, sigma=2), from=0, to=15,
   col="red", ylab="Hazard function", las=1)


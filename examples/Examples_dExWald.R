## The probability density function

curve(dExWALD(x, mu=1, sigma=4, nu=1), from=0, to=10, 
      ylim=c(0, 0.3), col="red", las=1, ylab="f(x)")

## The cumulative distribution and the Reliability function
curve(pExWALD(x, mu=1, sigma=4, nu=1), from=0, to=10, 
      col="red", las=1, ylab="F(x)")

curve(pExWALD(x, mu=1, sigma=4, nu=1, lower.tail=FALSE), 
      from=0, to=10, col="red", las=1, ylab="R(x)")

## The quantile function
p <- seq(from=0.01, to=0.99, length.out=100)
plot(x=qExWALD(p, mu=1.5, sigma=1.5, nu=2), y=p, xlab="Quantile",
     las=1, ylab="Probability")
curve(pExWALD(x, mu=1.5, sigma=1.5, nu=2), from=0, add=TRUE, col="red")

## The random function
x <- rExWALD(n=1000, mu=1.5, sigma=1.5, nu=2)
hist(x, freq=FALSE, xlab="x", 
     ylim=c(0, 0.35), las=1, main="")
curve(dExWALD(x, mu=1.5, sigma=1.5, nu=2), 
      from=0, to=20, add=TRUE, col="red")



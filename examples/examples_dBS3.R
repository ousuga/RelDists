# Example 1
# Plotting the mass function for different parameter values
curve(dBS3(x, mu=2, sigma=0.2), 
      from=0.001, to=10,
      ylim=c(0, 0.4), 
      col="royalblue1", lwd=2, 
      main="Density function",
      xlab="x", ylab="f(x)")
curve(dBS3(x, mu=2, sigma=0.4),
      col="tomato", 
      lwd=2,
      add=TRUE)
legend("topright", legend=c("mu=2, sigma=0.2", 
                            "mu=2, sigma=0.4"),
       col=c("royalblue1", "tomato"), lwd=2, cex=0.6)

# Example 2
# Checking if the cumulative curves converge to 1
curve(pBS3(x, mu=2, sigma=0.2), 
      from=0.00001, to=10,
      ylim=c(0, 1), 
      col="royalblue1", lwd=2, 
      main="Cumulative Distribution Function",
      xlab="x", ylab="f(x)")
curve(pBS3(x, mu=2, sigma=0.4),
      col="tomato", 
      lwd=2,
      add=TRUE)
legend("bottomright", legend=c("mu=2, sigma=0.2", 
                               "mu=2, sigma=0.4"),
       col=c("royalblue1", "tomato", "seagreen"), lwd=2, cex=0.5)

# Example 3
# The quantile function
p <- seq(from=0, to=0.999, length.out=100)
plot(x=qBS3(p, mu=2, sigma=0.2), y=p, xlab="Quantile",
     las=1, ylab="Probability", main="Quantile function ")
curve(pBS3(x, mu=2, sigma=0.2), 
      from=0, add=TRUE, col="tomato", lwd=2.5)

# Example 4
# The random function
x <- rBS3(n=10000, mu=2, sigma=0.2)
hist(x, freq=FALSE)
curve(dBS3(x, mu=2, sigma=0.2), from=0, to=10, 
      add=TRUE, col="tomato", lwd=2)

# Example 5
# The Hazard function
curve(hBS3(x, mu=2, sigma=0.2), from=0.001, to=4,
      col="tomato", ylab="Hazard function", las=1)


#Example 1
#Plotting the mass function for different parameter values
curve(dBS2(x, mu=1.0, sigma=100), 
      from=0.001, to=5,
      ylim=c(0, 3), 
      col="royalblue1", lwd=2, 
      main="Density function",
      xlab="x", ylab="f(x)")
curve(dBS2(x, mu=1.5, sigma=100),
      col="tomato", 
      lwd=2,
      add=TRUE)
curve(dBS2(x, mu=2.0, sigma=100),
      col="seagreen",
      lwd=2,
      add=TRUE)
legend("topright", legend=c("mu=1.0, sigma=100", 
                            "mu=1.5, sigma=100",
                            "mu=2.0, sigma=100"),
       col=c("royalblue1", "tomato", "seagreen"), lwd=2, cex=0.6)


curve(dBS2(x, mu=1, sigma=2), 
      from=0.001, to=2,
      ylim=c(0, 1.1), 
      col="royalblue1", lwd=2, 
      main="Density function",
      xlab="x", ylab="f(x)")
curve(dBS2(x, mu=1, sigma=5),
      col="tomato", 
      lwd=2,
      add=TRUE)
curve(dBS2(x, mu=1, sigma=10),
      col="seagreen",
      lwd=2,
      add=TRUE)
legend("topright", legend=c("mu=1, sigma=2", 
                            "mu=1, sigma=5",
                            "mu=1, sigma=10"),
       col=c("royalblue1", "tomato", "seagreen"), lwd=2, cex=0.6)


# Example 2
# Checking if the cumulative curves converge to 1
curve(pBS2(x, mu=0.5, sigma=0.5), 
      from=0.001, to=15,
      ylim=c(0, 1), 
      col="royalblue1", lwd=2, 
      main="Cumulative Distribution Function",
      xlab="x", ylab="f(x)")
curve(pBS2(x, mu=1, sigma=0.5),
      col="tomato", 
      lwd=2,
      add=TRUE)
curve(pBS2(x, mu=1.5, sigma=0.5),
      col="seagreen",
      lwd=2,
      add=TRUE)
legend("bottomright", legend=c("mu=0.5, sigma=0.5", 
                               "mu=1.0, sigma=0.5",
                               "mu=1.5, sigma=0.5"),
       col=c("royalblue1", "tomato", "seagreen"), lwd=2, cex=0.5)

# Example 3
# The quantile function
p <- seq(from=0, to=0.999, length.out=100)
plot(x=qBS2(p, mu=2.3, sigma=1.7), y=p, xlab="Quantile",
     las=1, ylab="Probability", main="Quantile function ")
curve(pBS2(x, mu=2.3, sigma=1.7), 
      from=0, add=TRUE, col="tomato", lwd=2.5)

# Example 4
# The random function
x <- rBS2(n=10000, mu=2.5, sigma=100)
hist(x, freq=FALSE)
curve(dBS2(x, mu=2.5, sigma=100), from=0, to=10, 
      add=TRUE, col="tomato", lwd=2)

# Example 5
# The Hazard function
curve(hBS2(x, mu=20, sigma=0.5), from=0.001, to=100,
      col="tomato", ylab="Hazard function", las=1)


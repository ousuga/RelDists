# Example 1
# Plotting the mass function for different parameter values
curve(dEXL(x, mu=0.5, sigma=0.5), 
      from=0, to=5,
      ylim=c(0, 1), 
      col="royalblue1", lwd=2, 
      main="Density function",
      xlab="x", ylab="f(x)")
curve(dEXL(x, mu=1, sigma=0.5),
      col="tomato", 
      lwd=2,
      add=TRUE)
curve(dEXL(x, mu=1.5, sigma=0.5),
      col="seagreen",
      lwd=2,
      add=TRUE)
legend("topright", legend=c("mu=0.5, sigma=0.5", 
                            "mu=1.0, sigma=0.5",
                            "mu=1.5, sigma=0.5"),
       col=c("royalblue1", "tomato", "seagreen"), lwd=2, cex=0.6)


curve(dEXL(x, mu=0.5, sigma=1), 
      from=0, to=5,
      ylim=c(0, 1), 
      col="royalblue1", lwd=2, 
      main="Density function",
      xlab="x", ylab="f(x)")
curve(dEXL(x, mu=1, sigma=1),
      col="tomato", 
      lwd=2,
      add=TRUE)
curve(dEXL(x, mu=1.5, sigma=1),
      col="seagreen",
      lwd=2,
      add=TRUE)
legend("topright", legend=c("mu=0.5, sigma=1", 
                            "mu=1.0, sigma=1",
                            "mu=1.5, sigma=1"),
       col=c("royalblue1", "tomato", "seagreen"), lwd=2, cex=0.6)


curve(dEXL(x, mu=0.5, sigma=1.5), 
      from=0., to=8,
      ylim=c(0, 1), 
      col="royalblue1", lwd=2, 
      main="Density function",
      xlab="x", ylab="f(x)")
curve(dEXL(x, mu=1, sigma=1.5),
      col="tomato", 
      lwd=2,
      add=TRUE)
curve(dEXL(x, mu=1.5, sigma=1.5),
      col="seagreen",
      lwd=2,
      add=TRUE)
legend("topright", legend=c("mu=0.5, sigma=1.5", 
                            "mu=1.0, sigma=1.5",
                            "mu=1.5, sigma=1.5"),
       col=c("royalblue1", "tomato", "seagreen"), lwd=2, cex=0.6)

# Example 2
# Checking if the cumulative curves converge to 1
curve(pEXL(x, mu=0.5, sigma=0.5), 
      from=0, to=5,
      ylim=c(0, 1), 
      col="royalblue1", lwd=2, 
      main="Cumulative Distribution Function",
      xlab="x", ylab="f(x)")
curve(pEXL(x, mu=1, sigma=0.5),
      col="tomato", 
      lwd=2,
      add=TRUE)
curve(pEXL(x, mu=1.5, sigma=0.5),
      col="seagreen",
      lwd=2,
      add=TRUE)
legend("bottomright", legend=c("mu=0.5, sigma=0.5", 
                               "mu=1.0, sigma=0.5",
                               "mu=1.5, sigma=0.5"),
       col=c("royalblue1", "tomato", "seagreen"), lwd=2, cex=0.5)


# Example 3
# The quantile function
p <- seq(from=0, to=0.99999, length.out=100)
plot(x=qEXL(p, mu=2.3, sigma=1.7), y=p, xlab="Quantile",
     las=1, ylab="Probability", main="Quantile function ")
curve(pEXL(x, mu=2.3, sigma=1.7), 
      from=0, add=TRUE, col="tomato", lwd=2.5)

# Comparing quantile, density and cumulative
p <- c(0.25, 0.5, 0.75)
quantile <- qEXL(p=p, mu=2.3, sigma=1.7) 

for(i in quantile){
  print(integrate(dEXL, lower=0, upper=i, mu=2.3, sigma=1.7))
}

pEXL(q=quantile, mu=2.3, sigma=1.7)

# Example 4
# The random function
x <- rEXL(n=10000, mu=1.5, sigma=2.5)
hist(x, freq=FALSE)
curve(dEXL(x, mu=1.5, sigma=2.5), from=0, to=20, 
      add=TRUE, col="tomato", lwd=2)

# Example 5
# The Hazard function
curve(hEXL(x, mu=1.5, sigma=2), from=0.001, to=4,
      col="tomato", ylab="Hazard function", las=1)

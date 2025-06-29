#Example 1
#Plotting the mass function for different parameter values
curve(dBS(x, mu=0.5, sigma=0.5), 
      from=0.001, to=5,
      ylim=c(0, 2), 
      col="royalblue1", lwd=2, 
      main="Density function",
      xlab="x", ylab="f(x)")
curve(dBS(x, mu=1, sigma=0.5),
      col="tomato", 
      lwd=2,
      add=TRUE)
curve(dBS(x, mu=1.5, sigma=0.5),
      col="seagreen",
      lwd=2,
      add=TRUE)
legend("topright", legend=c("mu=0.5, sigma=0.5", 
                            "mu=1.0, sigma=0.5",
                            "mu=1.5, sigma=0.5"),
       col=c("royalblue1", "tomato", "seagreen"), lwd=2, cex=0.6)


curve(dBS(x, mu=0.5, sigma=1), 
      from=0.001, to=5,
      ylim=c(0, 1.5), 
      col="royalblue1", lwd=2, 
      main="Density function",
      xlab="x", ylab="f(x)")
curve(dBS(x, mu=1, sigma=1),
      col="tomato", 
      lwd=2,
      add=TRUE)
curve(dBS(x, mu=1.5, sigma=1),
      col="seagreen",
      lwd=2,
      add=TRUE)
legend("topright", legend=c("mu=0.5, sigma=1", 
                            "mu=1.0, sigma=1",
                            "mu=1.5, sigma=1"),
       col=c("royalblue1", "tomato", "seagreen"), lwd=2, cex=0.6)


curve(dBS(x, mu=0.5, sigma=1.5), 
      from=0.001, to=8,
      ylim=c(0, 2), 
      col="royalblue1", lwd=2, 
      main="Density function",
      xlab="x", ylab="f(x)")
curve(dBS(x, mu=1, sigma=1.5),
      col="tomato", 
      lwd=2,
      add=TRUE)
curve(dBS(x, mu=1.5, sigma=1.5),
      col="seagreen",
      lwd=2,
      add=TRUE)
legend("topright", legend=c("mu=0.5, sigma=1.5", 
                            "mu=1.0, sigma=1.5",
                            "mu=1.5, sigma=1.5"),
       col=c("royalblue1", "tomato", "seagreen"), lwd=2, cex=0.6)

# Example 2
# Checking if the cumulative curves converge to 1
curve(pBS(x, mu=0.5, sigma=0.5), 
      from=0.001, to=5,
      ylim=c(0, 1), 
      col="royalblue1", lwd=2, 
      main="Cumulative Distribution Function",
      xlab="x", ylab="f(x)")
curve(pBS(x, mu=1, sigma=0.5),
      col="tomato", 
      lwd=2,
      add=TRUE)
curve(pBS(x, mu=1.5, sigma=0.5),
      col="seagreen",
      lwd=2,
      add=TRUE)
legend("bottomright", legend=c("mu=0.5, sigma=0.5", 
                               "mu=1.0, sigma=0.5",
                               "mu=1.5, sigma=0.5"),
       col=c("royalblue1", "tomato", "seagreen"), lwd=2, cex=0.5)
curve(pBS(x, mu=0.5, sigma=0.5, lower.tail=FALSE), 
      from=0.001, to=5,
      ylim=c(0, 1), 
      col="royalblue1", lwd=2, 
      main="Cumulative Distribution Function",
      xlab="x", ylab="f(x)")
curve(pBS(x, mu=1, sigma=0.5, lower.tail=FALSE),
      col="tomato", 
      lwd=2,
      add=TRUE)
curve(pBS(x, mu=1.5, sigma=0.5, lower.tail=FALSE),
      col="seagreen",
      lwd=2,
      add=TRUE)
legend("topright", legend=c("mu=0.5, sigma=0.5", 
                            "mu=1.0, sigma=0.5",
                            "mu=1.5, sigma=0.5"),
       col=c("royalblue1", "tomato", "seagreen"), lwd=2, cex=0.5)

#example 3
## The quantile function
p <- seq(from=0, to=0.999, length.out=100)
plot(x=qBS(p, mu=2.3, sigma=1.7), y=p, xlab="Quantile",
     las=1, ylab="Probability", main="Quantile function ")
curve(pBS(x, mu=2.3, sigma=1.7), 
      from=0, add=TRUE, col="tomato", lwd=2.5)

#some values
p <- c(0.25, 0.5, 0.75)
quantile <- qBS(p=p, mu=2.3, sigma=1.7) 
for(i in quantile){
  print(integrate(dBS, lower=0, upper=i, mu=2.3, sigma=1.7))
}

#example 4
## The random function
x <- rBS(n=10000, mu=20, sigma=0.5)
hist(x, freq=FALSE)
curve(dBS(x, mu=20, sigma=0.5), from=0, to=100, 
      add=TRUE, col="tomato", lwd=2)

#example 5
## The Hazard function
curve(hBS(x, mu=20, sigma=0.5), from=0.001, to=100,
      col="tomato", ylab="Hazard function", las=1)


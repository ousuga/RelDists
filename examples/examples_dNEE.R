# Example 1
# Plotting the mass function for different parameter values
curve(dNEE(x, mu=0.2, sigma=0.3),
      from=0, to=8, col="cadetblue3", las=1, ylab="f(x)")

curve(dNEE(x, mu=1, sigma=4),
      add=TRUE, col= "purple")

curve(dNEE(x, mu=1.5, sigma=22),
      add=TRUE, col="goldenrod")

curve(dNEE(x, mu=0.5, sigma=2),
      add=TRUE, col="green3")

legend("topright", col=c("cadetblue3", "purple", "goldenrod", "green3"), lty=1, bty="n",
       legend=c("mu=0.2, sigma=0.3",
                "mu=1.0, sigma=4",
                "mu=1.5, sigma=22",
                "mu=0.5, sigma=2"))

# Example 2
# Checking if the cumulative curves converge to 1
curve(pNEE(x, mu=0.2, sigma=0.3), ylim=c(0, 1),
      from=0, to=8, col="cadetblue3", las=1, ylab="F(x)")

curve(pNEE(x, mu=1, sigma=4),
      add=TRUE, col= "purple")

curve(pNEE(x, mu=1.5, sigma=22),
      add=TRUE, col="goldenrod")

curve(pNEE(x, mu=0.5, sigma=2),
      add=TRUE, col="green3")

legend("bottomright", col=c("cadetblue3", "purple", "goldenrod", "green3"), lty=1, bty="n",
       legend=c("mu=0.2, sigma=0.3",
                "mu=1.0, sigma=4",
                "mu=1.5, sigma=22",
                "mu=0.5, sigma=2"))

# Example 3
# Checking the quantile function
mu <- 0.5
sigma <- 2
p <- seq(from=0, to=0.999, length.out=100)
plot(x=qNEE(p, mu=mu, sigma=sigma), y=p, xlab="Quantile",
     las=1, ylab="Probability")
curve(pNEE(x, mu=mu, sigma=sigma), from=0, add=TRUE, col="red")

# Example 4
# Comparing the random generator output with
# the theoretical probabilities
mu <- 0.5
sigma <- 2
x <- rNEE(n=10000, mu=mu, sigma=sigma)
hist(x, freq=FALSE)
curve(dNEE(x, mu=mu, sigma=sigma), col="tomato", add=TRUE)





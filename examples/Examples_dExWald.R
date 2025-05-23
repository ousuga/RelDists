# Example 1
# Plotting the mass function for different parameter values
curve(dExWALD(x, mu=0.15, sigma=52.5, nu=50), ylim=c(0, 0.005),
      from=0, to=1200, col="cadetblue3", las=1, ylab="f(x)")

curve(dExWALD(x, mu=0.20, sigma=70, nu=50),
      add=TRUE, col= "purple")

curve(dExWALD(x, mu=0.25, sigma=87.5, nu=50),
      add=TRUE, col="goldenrod")

curve(dExWALD(x, mu=0.20, sigma=70, nu=115),
      add=TRUE, col="tomato")

curve(dExWALD(x, mu=0.20, sigma=70, nu=35),
      add=TRUE, col="blue")

legend("topright", col=c("cadetblue3", "purple", "goldenrod",
                         "tomato", "blue"), 
       lty=1, bty="n",
       legend=c("mu=0.15, sigma=52.5, nu=50",
                "mu=0.20, sigma=70.0, nu=50",
                "mu=0.25, sigma=87.5, nu=50",
                "mu=0.20, sigma=70.0, nu=115",
                "mu=0.20, sigma=70.0, nu=35"))

# Example 2
# Checking if the cumulative curves converge to 1
curve(pExWALD(x, mu=0.15, sigma=52.5, nu=50), ylim=c(0, 1),
      from=0, to=1200, col="cadetblue3", las=1, ylab="F(x)")

curve(pExWALD(x, mu=0.20, sigma=70, nu=50),
      add=TRUE, col= "purple")

curve(pExWALD(x, mu=0.25, sigma=87.5, nu=50),
      add=TRUE, col="goldenrod")

curve(pExWALD(x, mu=0.20, sigma=70, nu=115),
      add=TRUE, col="tomato")

curve(pExWALD(x, mu=0.20, sigma=70, nu=35),
      add=TRUE, col="blue")

legend("bottomright", col=c("cadetblue3", "purple", "goldenrod",
                         "tomato", "blue"), 
       lty=1, bty="n",
       legend=c("mu=0.15, sigma=52.5, nu=50",
                "mu=0.20, sigma=70.0, nu=50",
                "mu=0.25, sigma=87.5, nu=50",
                "mu=0.20, sigma=70.0, nu=115",
                "mu=0.20, sigma=70.0, nu=35"))

# Example 3
# Checking the quantile function
mu <- 5
sigma <- 3
nu <- 2
p <- seq(from=0.1, to=0.99, length.out=100)
plot(x=qExWALD(p, mu=mu, sigma=sigma, nu=nu), y=p, xlab="Quantile",
     las=1, ylab="Probability")
curve(pExWALD(x, mu=mu, sigma=sigma, nu=nu), from=0, add=TRUE, col="red")

# Example 4
# Comparing the random generator output with
# the theoretical probabilities
mu <- 0.2
sigma <- 70
nu <- 35
x <- rExWALD(n=10000, mu=mu, sigma=sigma, nu=nu)
hist(x, freq=FALSE)
curve(dExWALD(x, mu=mu, sigma=sigma, nu=nu), col="tomato", add=TRUE)


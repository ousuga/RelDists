## The probability density function 
curve(dWGEE(x, mu = 5, sigma = 0.5, nu = 1), from = 0, to = 6, 
ylim = c(0, 1), col = "red", las = 1, ylab = "The probability density function")

## The cumulative distribution and the Reliability function
par(mfrow = c(1, 2))
curve(pWGEE(x, mu = 5, sigma = 0.5, nu = 1), from = 0, to = 6, 
ylim = c(0, 1), col = "red", las = 1, ylab = "The cumulative distribution function")
curve(pWGEE(x, mu = 5, sigma = 0.5, nu = 1, lower.tail = FALSE), 
from = 0, to = 6, ylim = c(0, 1), col = "red", las = 1, ylab = "The Reliability function")

## The quantile function
p <- seq(from = 0, to = 0.99999, length.out = 100)
plot(x = qWGEE(p = p, mu = 5, sigma = 0.5, nu = 1), y = p, 
xlab = "Quantile", las = 1, ylab = "Probability")
curve(pWGEE(x, mu = 5, sigma = 0.5, nu = 1), from = 0, add = TRUE, 
col = "red")

## The random function
hist(rWGEE(1000, mu = 5, sigma = 0.5, nu = 1), freq = FALSE, xlab = "x", 
ylim = c(0, 1), las = 1, main = "")
curve(dWGEE(x, mu = 5, sigma = 0.5, nu = 1),  from = 0, add = TRUE, 
col = "red", ylim = c(0, 1))

## The Hazard function(
par(mfrow=c(1,1))
curve(hWGEE(x, mu = 5, sigma = 0.5, nu = 1), from = 0, to = 6, 
ylim = c(0, 1.4), col = "red", ylab = "The hazard function", las = 1)

## The probability density function 
curve(dWG(x, mu = 0.9, sigma = 2, nu = 0.5), from = 0, to = 3, 
ylim = c(0, 1.1), col = "red", las = 1, ylab = "f(x)")

## The cumulative distribution and the Reliability function
par(mfrow = c(1, 2))
curve(pWG(x, mu = 0.9, sigma = 2, nu = 0.5), from = 0, to = 3, 
ylim = c(0, 1), col = "red", las = 1, ylab = "F(x)")
curve(pWG(x, mu = 0.9, sigma = 2, nu = 0.5, lower.tail = FALSE), 
from = 0, to = 3, ylim = c(0, 1), col = "red", las = 1, ylab = "R(x)")

## The quantile function
p <- seq(from = 0, to = 0.99999, length.out = 100)
plot(x = qWG(p = p, mu = 0.9, sigma = 2, nu = 0.5), y = p, 
xlab = "Quantile", las = 1, ylab = "Probability")
curve(pWG(x,mu = 0.9, sigma = 2, nu = 0.5), from = 0, add = TRUE, 
col = "red")

## The random function
hist(rWG(1000, mu = 0.9, sigma = 2, nu = 0.5), freq = FALSE, xlab = "x", 
ylim = c(0, 1.8), las = 1, main = "")
curve(dWG(x, mu = 0.9, sigma = 2, nu = 0.5),  from = 0, add = TRUE, 
col = "red", ylim = c(0, 1.8))

## The Hazard function(
par(mfrow=c(1,1))
curve(hWG(x, mu = 0.9, sigma = 2, nu = 0.5), from = 0, to = 8, 
ylim = c(0, 12), col = "red", ylab = "Hazard function", las = 1)


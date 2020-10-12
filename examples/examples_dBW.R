## The probability density function 
curve(dBW(x, mu = (1/4), sigma =1, nu=1, tau=2), from = 0, to = 10, 
      col = "red", las = 1, ylab = "f(x)")

## The cumulative distribution and the Reliability function
par(mfrow = c(1, 2))
curve(pBW(x, mu = (1/4), sigma =1, nu=1, tau=2), from = 0, to = 6, 
ylim = c(0, 1), col = "red", las = 1, ylab = "F(x)")
curve(pBW(x, mu = (1/4), sigma =1, nu=1, tau=2, lower.tail = FALSE), 
from = 0, to = 6, ylim = c(0, 1), col = "red", las = 1, ylab = "R(x)")

## The quantile function
p <- seq(from = 0, to = 0.99999, length.out = 100)
plot(x = qBW(p = p, mu = (1/4), sigma =1, nu=1, tau=2), y = p, 
xlab = "Quantile", las = 1, ylab = "Probability")
curve(pBW(x, mu = (1/4), sigma =1, nu=1, tau=2), from = 0, add = TRUE, 
col = "red")

## The random function
hist(rBW(1000, mu = (1/4), sigma =1, nu=1, tau=2), freq = FALSE, xlab = "x", 
ylim = c(0, 0.5), las = 1, main = "")
curve(dBW(x, mu = (1/4), sigma =1, nu=1, tau=2),  from = 0, add = TRUE, 
col = "red", ylim = c(0, 0.5))

## The Hazard function(
par(mfrow=c(1,1))
curve(hBW(x, mu = (1/4), sigma =1, nu=1, tau=2), from = 0, to = 2, 
ylim = c(0, 1), col = "red", ylab = "Hazard function", las = 1)

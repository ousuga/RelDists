## The probability density function 
curve(dGammaW(x, mu = 0.5, sigma = 2, nu=1), from = 0, to = 6, 
      col = "red", las = 1, ylab = "f(x)")

## The cumulative distribution and the Reliability function
par(mfrow = c(1, 2))
curve(pGammaW(x, mu = 0.5, sigma = 2, nu=1), from = 0, to = 3, 
ylim = c(0, 1), col = "red", las = 1, ylab = "F(x)")
curve(pGammaW(x, mu = 0.5, sigma = 2, nu=1, lower.tail = FALSE), 
from = 0, to = 3, ylim = c(0, 1), col = "red", las = 1, ylab = "R(x)")

## The quantile function
p <- seq(from = 0, to = 0.99999, length.out = 100)
plot(x = qGammaW(p = p, mu = 0.5, sigma = 2, nu=1), y = p, 
xlab = "Quantile", las = 1, ylab = "Probability")
curve(pGammaW(x, mu = 0.5, sigma = 2, nu=1), from = 0, add = TRUE, 
col = "red")

## The random function
hist(rGammaW(1000, mu = 0.5, sigma = 2, nu=1), freq = FALSE, xlab = "x", 
ylim = c(0, 1), las = 1, main = "")
curve(dGammaW(x, mu = 0.5, sigma = 2, nu=1),  from = 0, add = TRUE, 
col = "red", ylim = c(0, 1))

## The Hazard function
par(mfrow=c(1,1))
curve(hGammaW(x, mu = 0.5, sigma = 2, nu=1), from = 0, to = 2, 
ylim = c(0, 1), col = "red", ylab = "Hazard function", las = 1)

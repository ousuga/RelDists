## The probability density function
curve(dSZMW(x, mu = 2, sigma = 1.5, nu = 0.2), from = 0, to = 2, 
      ylim = c(0, 1.7), col = "red", las = 1, ylab = "f(x)")
## The cumulative distribution and the Reliability function
par(mfrow = c(1, 2))
curve(pSZMW(x, mu = 2, sigma = 1.5, nu = 0.2), from = 0, to = 2, ylim = c(0, 1),
      col = "red", las = 1, ylab = "F(x)")
curve(pSZMW(x, mu = 2, sigma = 1.5, nu = 0.2, lower.tail = FALSE), from = 0,
      to = 2, ylim = c(0, 1), col = "red", las = 1, ylab = "R(x)")

## The quantile function
p <- seq(from = 0, to = 0.99999, length.out = 100)
plot(x = qSZMW(p = p, mu = 2, sigma = 1.5, nu = 0.2), y = p, xlab = "Quantile",
     las = 1, ylab = "Probability")
curve(pSZMW(x, mu = 2, sigma = 1.5, nu = 0.2), from = 0, add = TRUE, col = "red")

## The random function
hist(rSZMW(n = 1000, mu = 2, sigma = 1.5, nu = 0.2), freq = FALSE, xlab = "x",
     las = 1, main = "")
curve(dSZMW(x, mu = 2, sigma = 1.5, nu = 0.2),  from = 0, add = TRUE, col = "red")

## The Hazard function
par(mfrow=c(1,1))
curve(hSZMW(x, mu = 2, sigma = 1.5, nu = 0.2), from = 0, to = 3, ylim = c(0, 8),
      col = "red", ylab = "Hazard function", las = 1)

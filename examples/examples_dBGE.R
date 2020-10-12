## The probability density function 
curve(dBGE(x, mu = 1.5, sigma =1.7, nu=1, tau=1), from = 0, to = 3, 
      col = "red", las = 1, ylab = "f(x)")

## The cumulative distribution and the Reliability function
par(mfrow = c(1, 2))
curve(pBGE(x, mu = 1.5, sigma =1.7, nu=1, tau=1), from = 0, to = 6, 
      ylim = c(0, 1), col = "red", las = 1, ylab = "F(x)")
curve(pBGE(x, mu = 1.5, sigma =1.7, nu=1, tau=1, lower.tail = FALSE), 
      from = 0, to = 6, ylim = c(0, 1), col = "red", las = 1, ylab = "R(x)")

## The quantile function
p <- seq(from = 0, to = 0.99999, length.out = 100)
plot(x = qBGE(p = p, mu = 1.5, sigma =1.7, nu=1, tau=1), y = p, 
     xlab = "Quantile", las = 1, ylab = "Probability")
curve(pBGE(x, mu = (1/4), sigma =1, nu=1, tau=2), from = 0, add = TRUE, 
      col = "red")

## The random function
hist(rBGE(1000, mu = 1.5, sigma =1.7, nu=1, tau=1), freq = FALSE, xlab = "x", 
     ylim = c(0, 1), las = 1, main = "")
curve(dBGE(x, mu = 1.5, sigma =1.7, nu=1, tau=1),  from = 0, add = TRUE, 
      col = "red", ylim = c(0, 0.5))

## The Hazard function(
par(mfrow=c(1,1))
curve(hBGE(x, mu = 0.9, sigma =0.5, nu=1, tau=1), from = 0, to = 2, 
      col = "red", ylab = "Hazard function", las = 1)

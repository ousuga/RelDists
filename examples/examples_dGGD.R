old_par <- par(mfrow = c(1, 1)) # save previous graphical parameters

## The probability density function 
par(mfrow = c(1, 1))
curve(dGGD(x, mu=1, sigma=0.3, nu=1.5), from = 0, to = 4, 
      col = "red", las = 1, ylab = "f(x)")

## The cumulative distribution and the Reliability function
par(mfrow = c(1, 2))
curve(pGGD(x, mu=1, sigma=0.3, nu=1.5), from = 0, to = 4, 
      ylim = c(0, 1), col = "red", las = 1, ylab = "F(x)")
curve(pGGD(x, mu=1, sigma=0.3, nu=1.5, lower.tail = FALSE), 
      from = 0, to = 4, ylim = c(0, 1), col = "red", las = 1, ylab = "R(x)")

## The quantile function
p <- seq(from = 0, to = 0.99999, length.out = 100)
plot(x = qGGD(p=p, mu=1, sigma=0.3, nu=1.5), y = p, 
     xlab = "Quantile", las = 1, ylab = "Probability")
curve(pGGD(x, mu=1, sigma=0.3, nu=1.5), from = 0, add = TRUE, 
      col = "red")

## The random function
hist(rGGD(1000, mu=1, sigma=0.3, nu=1.5), freq = FALSE, xlab = "x", 
     las = 1, ylim = c(0, 0.7), main = "")
curve(dGGD(x,mu=1, sigma=0.3, nu=1.5), from = 0, to =8, add = TRUE, 
      col = "red")

## The Hazard function
par(mfrow=c(1,1))
curve(hGGD(x, mu=1, sigma=0.3, nu=1.5), from = 0, to = 3, col = "red",
      ylab = "The hazard function", las = 1)

par(old_par) # restore previous graphical parameters
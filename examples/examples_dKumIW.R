old_par <- par(mfrow = c(1, 1)) # save previous graphical parameters

## The probability density function 
par(mfrow = c(1, 1))
curve(dKumIW(x, mu = 1.5, sigma=  1.5, nu = 1), from = 0, to = 8.5, 
      col = "red", las = 1, ylab = "f(x)")

## The cumulative distribution and the Reliability function
par(mfrow = c(1, 2))
curve(pKumIW(x, mu = 1.5, sigma=  1.5, nu = 1), from = 0, to = 8.5, 
      ylim = c(0, 1), col = "red", las = 1, ylab = "F(x)")
curve(pKumIW(x, mu = 1.5, sigma=  1.5, nu = 1, lower.tail = FALSE), 
      from = 0, to = 6, ylim = c(0, 1), col = "red", las = 1, ylab = "R(x)")

## The quantile function
p <- seq(from = 0, to = 0.99999, length.out = 100)
plot(x = qKumIW(p=p, mu = 1.5, sigma=  1.5, nu = 10), y = p, 
     xlab = "Quantile", las = 1, ylab = "Probability")
curve(pKumIW(x, mu = 1.5, sigma=  1.5, nu = 10), from = 0, add = TRUE, 
      col = "red")

## The random function
hist(rKumIW(1000, mu = 1.5, sigma=  1.5, nu = 5), freq = FALSE, xlab = "x", 
     las = 1, ylim = c(0, 1.5), main = "")
curve(dKumIW(x, mu = 1.5, sigma=  1.5, nu = 5), from = 0, to =8, add = TRUE, 
      col = "red")

## The Hazard function
par(mfrow=c(1,1))
curve(hKumIW(x, mu = 1.5, sigma=  1.5, nu = 1), from = 0, to = 3, 
      ylim = c(0, 0.7), col = "red", ylab = "Hazard function", las = 1)

par(old_par) # restore previous graphical parameters
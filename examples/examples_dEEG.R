old_par <- par(mfrow = c(1, 1)) # save previous graphical parameters

## The probability density function 
par(mfrow=c(1,1))
curve(dEEG(x, mu = 1, sigma =3), from = 0, to = 10, 
      col = "red", las = 1, ylab = "f(x)")

## The cumulative distribution and the Reliability function
par(mfrow = c(1, 2))
curve(pEEG(x, mu = 1, sigma =3), from = 0, to = 10, 
      ylim = c(0, 1), col = "red", las = 1, ylab = "F(x)")
curve(pEEG(x, mu = 1, sigma =3, lower.tail = FALSE), 
      from = 0, to = 6, ylim = c(0, 1), col = "red", las = 1, ylab = "R(x)")

## The quantile function
p <- seq(from = 0, to = 0.99999, length.out = 100)
plot(x = qEEG(p = p, mu = 1, sigma =0.5), y = p, 
     xlab = "Quantile", las = 1, ylab = "Probability")
curve(pEEG(x, mu = 1, sigma =0.5), from = 0, add = TRUE, 
      col = "red")

## The random function
hist(rEEG(1000, mu = 1, sigma =1), freq = FALSE, xlab = "x", 
     ylim = c(0, 0.9), las = 1, main = "")
curve(dEEG(x, mu = 1, sigma =1),  from = 0, add = TRUE, 
      col = "red", ylim = c(0, 0.8))

## The Hazard function
par(mfrow=c(1,1))
curve(hEEG(x, mu = 1, sigma =0.5), from = 0, to = 2, 
      col = "red", ylab = "Hazard function", las = 1)

par(old_par) # restore previous graphical parameters
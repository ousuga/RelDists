## The probability density function 
curve(dEMWEx(x, mu = 49.046, sigma =3.148, nu=0.00005, tau=0.1), from=0, to=100,
      col = "red", las = 1, ylab = "f(x)")

## The cumulative distribution and the Reliability function
par(mfrow = c(1, 2))
curve(pEMWEx(x, mu = (1/4), sigma =1, nu=1, tau=2), from = 0, to = 1, 
      ylim = c(0, 1), col = "red", las = 1, ylab = "F(x)")
curve(pEMWEx(x, mu = (1/4), sigma =1, nu=1, tau=2, lower.tail = FALSE), 
      from = 0, to = 1, ylim = c(0, 1), col = "red", las = 1, ylab = "R(x)")

## The quantile function
p <- seq(from = 0, to = 0.99999, length.out = 100)
plot(x = qEMWEx(p = p, mu = 49.046, sigma =3.148, nu=0.00005, tau=0.1), y = p, 
     xlab = "Quantile", las = 1, ylab = "Probability")
curve(pEMWEx(x, mu = 49.046, sigma =3.148, nu=0.00005, tau=0.1), from = 0, add = TRUE, 
      col = "red")

## The random function
hist(rEMWEx(1000, mu = (1/4), sigma =1, nu=1, tau=2), freq = FALSE, xlab = "x", 
     las = 1, main = "")
curve(dEMWEx(x, mu = (1/4), sigma =1, nu=1, tau=2),  from = 0, add = TRUE, 
      col = "red", ylim = c(0, 0.5))

## The Hazard function(
par(mfrow=c(1,1))
curve(hEMWEx(x, mu = 49.046, sigma =3.148, nu=0.00005, tau=0.1), from = 0, to = 80, 
      col = "red", ylab = "Hazard function", las = 1)

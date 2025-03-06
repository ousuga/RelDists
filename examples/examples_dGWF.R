old_par <- par(mfrow = c(1, 1)) # save previous graphical parameters

## The probability density function
curve(
    dGWF(x, mu = 5, sigma = 2, nu = -0.2),
    from = 0, to = 5, col = "red", las = 1, ylab = "f(x)"
)

## The cumulative distribution and the Reliability function
par(mfrow = c(1, 2))
curve(
    pGWF(x, mu = 5, sigma = 2, nu = -0.2),
    from = 0, to = 5, ylim = c(0, 1),
    col = "red", las = 1, ylab = "F(x)"
)
curve(
    pGWF(
        x, mu = 5, sigma = 2, nu = -0.2, 
        lower.tail = FALSE
    ),
    from = 0, to = 5, ylim = c(0, 1),
    col = "red", las = 1, ylab = "R(x)"
)

## The quantile function
p <- seq(from = 0, to = 0.999, length.out = 100)
plot(
    x = qGWF(p, mu = 5, sigma = 2, nu = -0.2),
    y = p, xlab = "Quantile", las = 1,
    ylab = "Probability"
)
curve(
    pGWF(x, mu = 5, sigma = 2, nu = -0.2),
    from = 0, add = TRUE, col = "red"
)

## The random function
hist(
    rGWF(n = 10000, mu = 5, sigma = 2, nu = -0.2),
    freq = FALSE, xlab = "x", las = 1, main = "", ylim = c(0, 2.0)
)
curve(dGWF(x, mu = 5, sigma = 2, nu = -0.2),
    from = 0, add = TRUE, col = "red"
)

## The Hazard function
par(mfrow = c(1, 1))
curve(
    hGWF(x, mu = 0.003, sigma = 5e-6, nu = 0.025),
    from = 0, to = 250, col = "red", 
    ylab = "Hazard function", las = 1
)

par(old_par) # restore previous graphical parameters

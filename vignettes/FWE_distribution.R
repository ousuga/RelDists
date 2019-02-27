## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "##"
)

## ----pdfFWE1, echo=FALSE, message=FALSE----------------------------------
require(RelDists)
curve(dFWE(x, mu=0.75, sigma=0.5), from=0, to=3, ylim=c(0, 1.7), 
      col="red", ylab="f(x)", las=1)
curve(dFWE(x, mu=2.00, sigma=3.0), col="blue", add=TRUE)
curve(dFWE(x, mu=0.75, sigma=1.3), col="green", add=TRUE)

cap1 <- expression(paste(mu, "=", 0.75, ", ", sigma, "=", 0.5))
cap2 <- expression(paste(mu, "=", 2.00, ", ", sigma, "=", 3.0))
cap3 <- expression(paste(mu, "=", 0.75, ", ", sigma, "=", 1.3))

legend('topright', legend=c(cap1, cap2, cap3), 
       col=c("red", "blue", "green"), lty=c(1, 1, 1), bty="n")

## ----cumulativeFWE2, echo=FALSE------------------------------------------
curve(pFWE(x, mu=0.75, sigma=0.5), from=0, to=3, col="red",
      ylab="F(x)", las=1)
curve(pFWE(x, mu=2, sigma=3), col="blue", add=TRUE)
curve(pFWE(x, mu=0.75, sigma=1.3), col="green", add=TRUE)

cap1 <- expression(paste(mu, "=", 0.75, ", ", sigma, "=", 0.5))
cap2 <- expression(paste(mu, "=", 2.00, ", ", sigma, "=", 3.0))
cap3 <- expression(paste(mu, "=", 0.75, ", ", sigma, "=", 1.3))

legend('bottomright', legend=c(cap1, cap2, cap3), 
       col=c("red", "blue", "green"), lty=c(1, 1, 1), bty="n")

## ----FWE3, echo=FALSE----------------------------------------------------
curve(hFWE(x, mu=0.75, sigma=0.5), from=0, to=3, ylim=c(0, 3), 
      col="red", ylab="h(x)", las=1)
curve(hFWE(x, mu=2, sigma=3), col="blue", add=T)
curve(hFWE(x, mu=0.75, sigma=1.3), col="green", add=T)

cap1 <- expression(paste(mu, "=", 0.75, ", ", sigma, "=", 0.5))
cap2 <- expression(paste(mu, "=", 2.00, ", ", sigma, "=", 3.0))
cap3 <- expression(paste(mu, "=", 0.75, ", ", sigma, "=", 1.3))

legend('bottomright', legend=c(cap1, cap2, cap3), 
       col=c("red", "blue", "green"), lty=c(1, 1, 1), bty="n")

## ------------------------------------------------------------------------
y <- c(2.160, 0.746, 0.402, 0.954, 0.491, 6.560, 4.992, 0.347,
       0.150, 0.358, 0.101, 1.359, 3.465, 1.060, 0.614, 1.921,
       4.082, 0.199, 0.605, 0.273, 0.070, 0.062, 5.320)

## ----message=FALSE-------------------------------------------------------
require(RelDists)
require(gamlss)
mod <- gamlss(y~1, sigma.fo=~1, family='FWE',
              control=gamlss.control(n.cyc=5000, trace=FALSE))

## ------------------------------------------------------------------------
exp(coef(mod, what='mu'))
exp(coef(mod, what='sigma'))

## ----FWE4----------------------------------------------------------------
hist(y, freq=FALSE, breaks=10, ylim=c(0, 1.7), las=1, main='')
curve(dFWE(x, mu=0.21, sigma=0.26), from=0.01, add=TRUE, lwd=2)

## ------------------------------------------------------------------------
y <- rFWE(n=100, mu=0.75, sigma=0.5)

## ---- message=FALSE------------------------------------------------------
mod <- gamlss(y~1, sigma.fo=~1, family='FWE',
              control=gamlss.control(n.cyc=250, trace=FALSE))

## ------------------------------------------------------------------------
exp(coef(mod, what='mu'))
exp(coef(mod, what='sigma'))

## ----FWE5----------------------------------------------------------------
hist(y, freq=FALSE, breaks=10, ylim=c(0, 1.5), las=1, main='')
curve(dFWE(x, mu=exp(coef(mod, what='mu')), sigma=exp(coef(mod, what='sigma'))), 
      from=0.01, add=TRUE, lwd=2)

## ------------------------------------------------------------------------
n <- 200
x1 <- runif(n)
x2 <- runif(n)
mu <- exp(1.21 - 3 * x1)
sigma <- exp(1.26 - 2 * x2)
y <- rFWE(n=n, mu, sigma)

## ------------------------------------------------------------------------
mod <- gamlss(y~x1, sigma.fo=~x2, family=FWE,
              control=gamlss.control(n.cyc=5000, trace=FALSE))

## ------------------------------------------------------------------------
coef(mod, what="mu")
coef(mod, what="sigma")


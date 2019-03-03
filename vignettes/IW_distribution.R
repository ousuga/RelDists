
  
## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "##"
)

## ---- IW1, echo=FALSE, message=FALSE-------------------------------------
library(RelDists)
curve(dIW(x, mu=2.0, sigma=1.0), from=0, to=10, las=1, xlab="Density", ylim=c(0,0.6), col="red", ylab="f(x)")
curve(dIW(x, mu=1.2, sigma=1.5), col= "purple", add=TRUE)
curve(dIW(x, mu=5.0, sigma=2.5), col= "green3", add=TRUE)

cap1 <- expression(paste(mu, "=", 2, ", ", sigma, "=", 1))
cap2 <- expression(paste(mu, "=", 1.2, ", ", sigma, "=", 1.5))
cap3 <- expression(paste(mu, "=", 5, ", ", sigma, "=", 2.5))

legend("topright", legend=c(cap1, cap2, cap3), 
       col=c("red","purple", "green3"), lty=c(1,1,1), bty="n")

## ---- IW2, echo=FALSE----------------------------------------------------

curve(pIW(x, mu=2.0, sigma=1.0), from=0, to=10, las=1, ylim=c(0,1), col="red", ylab="F(x)")
curve(pIW(x, mu=1.2, sigma=1.5), col= "purple", add=TRUE)
curve(pIW(x, mu=5.0, sigma=2.5), col= "green3", add=TRUE)

cap1 <- expression(paste(mu, "=", 2, ", ", sigma, "=", 1))
cap2 <- expression(paste(mu, "=", 1.2, ", ", sigma, "=", 1.5))
cap3 <- expression(paste(mu, "=", 5, ", ", sigma, "=", 2.5))

legend("bottomright", legend=c(cap1, cap2, cap3), 
       col=c("red","purple", "green3"), lty=c(1,1,1), bty="n")

## ----IW3, echo=FALSE-----------------------------------------------------
curve(hIW(x, mu=1, sigma=0.6), from=0, to=15, las = 1, ylim = c(0, 1), col="red", ylab="h(x)")
curve(hIW(x, mu=2, sigma=1), col="purple", add=TRUE)
curve(hIW(x, mu=1.2, sigma=1.5), col = "green3", add=TRUE)

cap1 <- expression(paste(mu, "=", 1, ", ", sigma, "=", 0.6))
cap2 <- expression(paste(mu, "=", 2, ", ", sigma, "=", 1))
cap3 <- expression(paste(mu, "=", 1.2, ", ", sigma, "=", 1.5))

legend("topright", legend=c(cap1, cap2, cap3), 
       col=c("red","purple", "green3"), lty=c(1,1,1), bty="n")

## ------------------------------------------------------------------------

y <- c(23, 32, 35, 20, 26, 24, 24, 14,
       13, 16, 5, 11, 5, 12, 12, 7,
       6, 6, 9, 9, 11, 12, 25, 26,
       15, 12, 12, 7, 15, 12, 29, 10,
       7, 10, 15, 20, 20, 17, 24, 31,
       26, 9, 16, 14, 18, 16, 14, 12
)

## ----message=FALSE-------------------------------------------------------
require(RelDists)
require(gamlss)
mod <- gamlss(y~1, mu.fo=~1, sigma.fo=~1, family='IW',
              control=gamlss.control(n.cyc=1500, trace=FALSE))

## ------------------------------------------------------------------------
exp(coef(mod, what='mu'))
exp(coef(mod, what='sigma'))

## ----IW4-----------------------------------------------------------------
hist(y, freq=FALSE, breaks=15, main='')
curve(dIW(x, mu=exp(coef(mod, what='mu')), 
          sigma=exp(coef(mod, what='sigma'))),
      from=0.01, add=TRUE, lwd=2)

## ------------------------------------------------------------------------

y <- rIW(n=100, mu=5, sigma=2.5)

## ---- message=FALSE------------------------------------------------------

mod <- gamlss(y~1, mu.fo=~1, sigma.fo=~1, family='IW',
              control=gamlss.control(n.cyc=250, trace=FALSE))

## ------------------------------------------------------------------------
exp(coef(mod, what='mu'))
exp(coef(mod, what='sigma'))

## ----IW5-----------------------------------------------------------------
hist(y, freq=FALSE, main='', ylim=c(0, 0.7))
curve(dIW(x, mu=exp(coef(mod, what='mu')), 
          sigma=exp(coef(mod, what='sigma'))), 
      from=0.01, add=TRUE, lwd=2)

## ----IW6-----------------------------------------------------------------
dt <- data.frame(y=c(12,15, 22, 24, 24, 32, 32, 33, 34, 38, 38, 43, 
                     44, 48, 52, 53, 54, 54, 55, 56, 57, 58, 58, 59,
                     60, 60, 60, 60, 61, 62, 63, 65, 65, 67, 68, 70,
                     70, 72, 73, 75, 76, 76, 81, 83, 84, 85, 87, 91,
                     95, 96, 98, 99, 109, 110, 121, 127,129, 131, 143, 146,
                     146, 175, 175, 211, 233, 258, 258, 263, 297, 341, 341, 376), 
                 
                 status=c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
                          TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
                          TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
                          TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
                          TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
                          FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE))

plot(density(dt$y), main='')
rug(dt$y)

## ------------------------------------------------------------------------
require(survival)
require(gamlss.cens)

mod <- gamlss(Surv(y, status) ~ 1,
              family=cens(IW), data=subset(dt, y < 40),
              control=gamlss.control(n.cyc=5000, trace=FALSE))

## ------------------------------------------------------------------------
exp(coef(mod, what='mu'))
exp(coef(mod, what='sigma'))

## ----IW7-----------------------------------------------------------------
curve(hIW(x, mu=exp(coef(mod, what='mu')), 
          sigma=exp(coef(mod, what='sigma'))),
      ylab='Hazard', xlab='Survival times (days)', from=0, to=500, las=1,add=TRUE, lwd=2)

## ------------------------------------------------------------------------
n <- 200
x1 <- rpois(n, lambda=2)
x2 <- runif(n)
mu <- exp(2 + -1 * x1)
sigma <- exp(2 - 2 * x2)
nu <- 2
y <- rIW(n=n, mu, sigma)

## ------------------------------------------------------------------------
mod <- gamlss(y~x1, sigma.fo=~x2, nu.fo=~1, family=IW,
              control=gamlss.control(n.cyc=5000, trace=FALSE))

## ------------------------------------------------------------------------
coef(mod, what="mu")
coef(mod, what="sigma")

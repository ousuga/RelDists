## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "##"
)

## ---- EW1, echo=FALSE, message=FALSE-------------------------------------
library(RelDists)
curve(dEW(x, mu=2.0, sigma=1.5, nu=0.5), from=0, to=2, las=1, xlab='Density',
      ylim=c(0,3), col="red", ylab="f(x)", las=1)
curve(dEW(x, mu=5.0, sigma=0.5, nu=10), col="green3", add=TRUE)
curve(dEW(x, mu=1.0, sigma=3.0, nu=1.0), col="purple", add=TRUE)

cap1 <- expression(paste(mu, "=", 2, ", ", sigma, "=", 1.5, ", ", nu, "=", 0.5))
cap2 <- expression(paste(mu, "=", 5, ", ", sigma, "=", 0.5, ", ", nu, "=", 10))
cap3 <- expression(paste(mu, "=", 1, ", ", sigma, "=", 3.0, ", ", nu, "=", 1.0))

legend('topright', legend=c(cap1, cap2, cap3), 
       col=c("red", "green3", "purple"), lty=c(1, 1, 1), bty="n")

## ---- EW2, echo=FALSE----------------------------------------------------
curve(pEW(x, mu=2.0, sigma=1.5, nu=0.5), 
      from=0, to=3,  col="red", ylab="F(x)", las=1)
curve(pEW(x, mu=5.0, sigma=0.5, nu=10), col="green", add=TRUE)
curve(pEW(x, mu=1.0, sigma=3.0, nu=1.0), col="purple", add=TRUE)  

cap1 <- expression(paste(mu, "=", 2, ", ", sigma, "=", 1.5, ", ", nu, "=", 0.5))
cap2 <- expression(paste(mu, "=", 5, ", ", sigma, "=", 0.5, ", ", nu, "=", 10))
cap3 <- expression(paste(mu, "=", 1, ", ", sigma, "=", 3.0, ", ", nu, "=", 1.0))

legend('bottomright', legend=c(cap1, cap2, cap3), 
       col=c("red", "green3", "purple"), lty=c(1, 1, 1), bty="n")

## ----EW3, echo=FALSE-----------------------------------------------------
curve(hEW(x, mu=2.0, sigma=1.5, nu=0.5), ylim=c(0, 7),
      from=0.001, to=3,  col="red", ylab="h(x)", las=1)
curve(hEW(x, mu=5.0, sigma=0.5, nu=10), col="green", add=TRUE)
curve(hEW(x, mu=1.0, sigma=3.0, nu=1.0), col="purple", add=TRUE)  

cap1 <- expression(paste(mu, "=", 2, ", ", sigma, "=", 1.5, ", ", nu, "=", 0.5))
cap2 <- expression(paste(mu, "=", 5, ", ", sigma, "=", 0.5, ", ", nu, "=", 10))
cap3 <- expression(paste(mu, "=", 1, ", ", sigma, "=", 3.0, ", ", nu, "=", 1.0))

legend('bottomright', legend=c(cap1, cap2, cap3), 
       col=c("red", "green3", "purple"), lty=c(1, 1, 1), bty="n")

## ------------------------------------------------------------------------
y <- c(1460, 4050, 3570, 2060, 1300, 1390, 1720, 6280, 1360, 7440,
       5320, 1400, 3240, 2710, 4520, 4840, 8320, 13900, 71500, 6250,
       2260, 318, 1330, 970, 1920, 15100, 2870, 20600, 3810, 726,
       7500, 7170, 2000, 829, 17300, 4740, 13400, 1940, 5660)

## ----message=FALSE-------------------------------------------------------
require(RelDists)
require(gamlss)
mod <- gamlss(y~1, sigma.fo=~1, nu.fo=~1, family='EW',
              control=gamlss.control(n.cyc=1500, trace=FALSE))

## ------------------------------------------------------------------------
exp(coef(mod, what='mu'))
exp(coef(mod, what='sigma'))
exp(coef(mod, what='nu'))

## ----EW4-----------------------------------------------------------------
hist(y, freq=FALSE, ylim=c(0, 0.0002), breaks=15, main='')
curve(dEW(x, mu=0.23, sigma=0.33, nu=19.03), from=0.01, add=TRUE, lwd=2)

## ------------------------------------------------------------------------
y <- rEW(n=100, mu=5, sigma=0.5, nu=10)

## ---- message=FALSE------------------------------------------------------
mod <- gamlss(y~1, sigma.fo=~1, nu.fo=~1, family='EW',
              control=gamlss.control(n.cyc=250, trace=FALSE))

## ------------------------------------------------------------------------
exp(coef(mod, what='mu'))
exp(coef(mod, what='sigma'))
exp(coef(mod, what='nu'))

## ----EW5-----------------------------------------------------------------
hist(y, freq=FALSE, main='', ylim=c(0, 1.8))
curve(dEW(x, mu=exp(coef(mod, what='mu')), 
          sigma=exp(coef(mod, what='sigma')), 
          nu=exp(coef(mod, what='nu'))), 
      from=0.01, add=TRUE, lwd=2)

## ----EW6-----------------------------------------------------------------
dt <- data.frame(y=c(7, 34, 42, 63, 64, 74, 83, 84, 91, 108, 112, 129, 
                     133, 133, 139, 140, 140, 146, 149, 154, 157, 160,
                     160, 165, 173, 176, 185, 218, 225, 241, 248, 273, 
                     277, 279, 297, 319, 405, 417, 420, 440, 523, 523, 
                     583, 594, 1101, 1116, 1146, 1226, 1349, 1412, 
                     1417) / 30.438, 
                 status=c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, 
                          TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, 
                          TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, 
                          TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, 
                          TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, 
                          FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, 
                          TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, 
                          FALSE, FALSE, TRUE))

plot(density(dt$y), main='')
rug(dt$y)

## ------------------------------------------------------------------------
require(survival)
require(gamlss.cens)

mod <- gamlss(Surv(y, status) ~ 1,
              family=cens(EW), data=subset(dt, y < 40),
              control=gamlss.control(n.cyc=5000, trace=FALSE))

## ------------------------------------------------------------------------
exp(coef(mod, what='mu'))
exp(coef(mod, what='sigma'))
exp(coef(mod, what='nu'))

## ----EW7-----------------------------------------------------------------
curve(hEW(x, mu=exp(coef(mod, what='mu')), 
          sigma=exp(coef(mod, what='sigma')), 
          nu=exp(coef(mod, what='nu'))),
      ylab='Hazard', xlab='Survival times (months)', from=0, to=50, las=1)

## ------------------------------------------------------------------------
n <- 200
x1 <- rpois(n, lambda=2)
x2 <- runif(n)
mu <- exp(2 + -3 * x1)
sigma <- exp(3 - 2 * x2)
nu <- 2
y <- rEW(n=n, mu, sigma, nu)

## ------------------------------------------------------------------------
mod <- gamlss(y~x1, sigma.fo=~x2, nu.fo=~1, family=EW,
              control=gamlss.control(n.cyc=5000, trace=FALSE))

## ------------------------------------------------------------------------
coef(mod, what="mu")
coef(mod, what="sigma")
coef(mod, what="nu")


#' @name OW
#' 
#' @title 
#' The Odd Weibull Distribution
#' 
#' @description 
#' Density, distribution function, quantile function, 
#' random generation and hazard function for the odd weibull distribution with
#' parameters \code{mu}, \code{theta} and \code{nu}.
#' 
#' @param x,q	vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. 
#' @param mu parameter one.    
#' @param theta parameter two.
#' @param nu parameter three.        
#' @param log,log.p	logical; if TRUE, probabilities p are given as log(p).	
#' @param lower.tail logical; if TRUE (default), probabilities are 
#' P[X <= x], otherwise, P[X > x].
#' 
#' @details 
#' The generalized power weibull with parameters \code{mu}, \code{theta}
#' and \code{nu} has density given by
#' 
#' f(x) = mu*theta*nu*x^(theta-1)*exp(mu*(x^theta))*(exp(mu*(x^theta))-1)^(nu-1)*(1+(exp(mu*(x^theta))-1)^nu)^-2
#'
#' for x > 0.
#' 
#' @return 
#' \code{dOW} gives the density, \code{pOW} gives the distribution 
#' function, \code{qOW} gives the quantile function, \code{rOW}
#' generates random deviates and \code{hOW} gives the hazard function.
#' 
#' @export
#' @examples  
#' ## The probability density function 
#' curve(dOW(x, mu = 2, theta = 3, nu = 0.2), from = 0, to = 4, ylim = c(0, 2), col = "red", las = 1, ylab = "The probability density function")
#' 
#' ## The cumulative distribution and the Reliability function
#' par(mfrow = c(1, 2))
#' curve(pOW(x, mu = 2, theta = 3, nu = 0.2), from = 0, to = 4, ylim = c(0, 1), col = "red", las = 1, ylab = "The cumulative distribution function")
#' curve(pOW(x, mu = 2, theta = 3, nu = 0.2, lower.tail = FALSE), from = 0, to = 4,  ylim = c(0, 1), col = "red", las = 1, ylab = "The Reliability function")
#' 
#' ## The quantile function
#' p <- seq(from = 0, to = 0.998, length.out = 100)
#' plot(x = qOW(p, mu = 2, theta = 3, nu = 0.2), y = p, xlab = "Quantile", las = 1, ylab = "Probability")
#' curve(pOW(x, mu = 2, theta = 3, nu = 0.2), from = 0, add = TRUE, col = "red")
#' 
#' ## The random function
#' hist(rOW(n = 10000, mu = 2, theta = 3, nu = 0.2), freq = FALSE, ylim = c(0, 2),xlab = "x", las = 1, main = "")
#' curve(dOW(x, mu = 2, theta = 3, nu = 0.2),  from = 0, ylim = c(0, 2), add = T, col = "red")
#' 
#' ## The Hazard function
#' curve(hOW(x, mu = 2, theta = 3, nu = 0.2), from = 0, to = 2.5, ylim = c(0, 3), col = "red", ylab = "The hazard function", las = 1)
OW <- function (mu.link = "log", sigma.link = "log", nu.link = "log") 
{
  mstats <- checklink("mu.link",    "Odd Weibull", substitute(mu.link),    c("log", "own"))
  dstats <- checklink("sigma.link", "Odd Weibull", substitute(sigma.link), c("log", "own"))
  vstats <- checklink("nu.link",    "Odd Weibull", substitute(nu.link),    c("log", "own"))
  
  structure(list(family = c("OW", "Odd Weibull"), 
                 parameters = list(mu = TRUE, sigma = TRUE, nu = TRUE), 
                 nopar = 3, 
                 type = "Continuous", 
                 
                 mu.link    = as.character(substitute(mu.link)), 
                 sigma.link = as.character(substitute(sigma.link)), 
                 nu.link    = as.character(substitute(nu.link)), 
                 
                 mu.linkfun    = mstats$linkfun, 
                 sigma.linkfun = dstats$linkfun, 
                 nu.linkfun    = vstats$linkfun, 
                 
                 mu.linkinv    = mstats$linkinv, 
                 sigma.linkinv = dstats$linkinv, 
                 nu.linkinv    = vstats$linkinv, 
                 
                 mu.dr    = mstats$mu.eta, 
                 sigma.dr = dstats$mu.eta, 
                 nu.dr    = vstats$mu.eta, 
                 
                 
                 dldm = function(y, mu, sigma, nu) {
                   exp1 <- exp(mu*y^sigma)
                   exp2 <- exp1-1
                   dexp1dm <- (y^sigma)*exp1                 
                   dldm <- 1/mu + y^sigma + ((nu-1)*dexp1dm)/exp2
                   dldm
                 },
                 
                 d2ldm2 = function(y, mu, sigma, nu) {
                   exp1 <- exp(mu*y^sigma)
                   exp2 <- exp1-1
                   dexp1dm <- (y^sigma)*exp1 
                   d2exp1dm2 <- y*dexp1dm
                   d2ldm2 <- -(-(1/mu^2) + ((nu-1)/exp2^2)*(exp2*d2exp1dm2-dexp1dm^2))^2
                   d2ldm2
                 }, 
                 
                 dldd = function(y, mu, sigma, nu) {
                   exp1 <- exp(mu*y^sigma)
                   exp2 <- exp1-1
                   dexp1dd <- mu*exp1*(y^sigma)*log(y)
                   dldd <- 1/sigma + log(y) + mu*(y^sigma)*log(y) + ((nu-1)*dexp1dd)/exp2
                   dldd
                 },
                 
                 d2ldd2 = function(y, mu, sigma, nu) {
                   exp1 <- exp(mu*y^sigma)
                   exp2 <- exp1-1
                   dexp1dm <- (y^sigma)*exp1
                   dexp1dd <- mu*exp1*(y^sigma)*log(y)
                   d2exp1dd2 <- mu*dexp1dm*(log(y))^2*(mu*(y^sigma)+1)
                   d2ldd2 <- -(-(1/sigma^2) + mu*(y^sigma)*(log(y))^2 + ((nu-1)/exp2^2)*(exp2*d2exp1dd2-dexp1dd^2))^2
                   d2ldd2
                 }, 
                 
                 dldv = function(y, mu, sigma, nu) {
                   exp1 <- exp(mu*y^sigma)
                   exp2 <- exp1-1
                   exp3 <- 1 + exp2^nu
                   dexp3dv <- (exp2^nu)*log(exp2)
                   dldv <- 1/nu + log(exp2) -(2*dexp3dv)/exp3
                   dldv 
                 }, 
                 
                 d2ldv2 = function(y, mu, sigma, nu) {
                   exp1 <- exp(mu*y^sigma)
                   exp2 <- exp1-1
                   exp3 <- 1 + exp2^nu
                   dexp3dv <- (exp2^nu)*log(exp2)
                   d2ldv2 <- -(-(1/nu^2) + (2*dexp3dv^2)/(exp3^2))^2
                   d2ldv2 
                 }, 
                 
                 d2ldmdd = function(y, mu, sigma, nu) {
                   exp1 <- exp(mu*y^sigma)
                   exp2 <- exp1-1
                   dexp1dm <- (y^sigma)*exp1
                   dexp1dd <- mu*exp1*(y^sigma)*log(y)
                   d2exp1dmdd <- (dexp1dd*(mu*(y^sigma)+1))/mu
                   d2ldmdd <- -((y^sigma)*log(y) + ((nu-1)/exp2^2)*(exp2*d2exp1dmdd - dexp1dm*dexp1dd))^2
                   d2ldmdd
                 },  
                 
                 d2ldmdv = function(y, mu, sigma, nu) {
                   exp1 <- exp(mu*y^sigma)
                   exp2 <- exp1-1
                   dexp1dm <- (y^sigma)*exp1
                   d2ldmdv <- -(dexp1dm/exp2)^2
                   d2ldmdv
                 },   
                 
                 d2ldddv = function(y, mu, sigma, nu) {
                   exp1 <- exp(mu*y^sigma)
                   exp2 <- exp1-1
                   dexp1dd <- mu*exp1*(y^sigma)*log(y)
                   d2ldddv <- -(dexp1dd/exp2)^2
                   d2ldddv 
                 },   
                 
                 
                 G.dev.incr = function(y, mu, sigma, nu, ...) -2*dOW(y, mu, sigma, nu, log = TRUE), 
                 rqres = expression(rqres(pfun = "pOW", type = "Continuous",  y = y, mu = mu, sigma = sigma, nu = nu)), 
                 
                 mu.initial = expression( mu <-  rep(0.5, length(y)) ), 
                 sigma.initial = expression( sigma <- rep(0.5, length(y)) ), 
                 nu.initial = expression( nu <- rep(0.5, length(y)) ), 
                 
                 mu.valid = function(mu) all(mu >  0), 
                 sigma.valid = function(sigma) all(sigma >  0), 
                 nu.valid = function(nu) all(nu > 0), 
                 
                 y.valid = function(y) all(y > 0)
  ), 
  class = c("gamlss.family", "family"))
}
#' @export
#' @rdname OW
dOW<-function(x,mu,theta,nu, log = FALSE){
  if (any(x<0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(mu<=0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(theta<=0)) 
    stop(paste("theta must be positive", "\n", ""))
  if (any(nu<=0)) 
    stop(paste("nu must be positive", "\n", ""))
  
  loglik<- log(mu) +log(theta) + log(nu) + (theta-1)*log(x) +
    mu*(x^theta) + (nu-1)*log(exp(mu*(x^theta))-1) -
    2*log(1+(exp(mu*(x^theta))-1)^nu)
  
  if (log == FALSE) 
    density<- exp(loglik)
  else 
    density <- loglik
  return(density)
}

#' @export
#' @rdname OW
pOW <- function(q,mu,theta,nu, lower.tail=TRUE, log.p = FALSE){
  if (any(q<0)) 
    stop(paste("q must be positive", "\n", ""))
  if (any(mu<=0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(theta<=0)) 
    stop(paste("theta must be positive", "\n", ""))
  if (any(nu<=0)) 
    stop(paste("nu must be positive", "\n", ""))
  
  cdf <- 1 - (1 + (exp(mu*(q^theta))-1)^nu )^(-1)
  if (lower.tail == TRUE) 
    cdf <- cdf
  else cdf <- 1 - cdf
  if (log.p == FALSE) 
    cdf <- cdf
  else cdf <- log(cdf)
  cdf
}

#' @export
#' @rdname OW
qOW <- function(p,mu,theta,nu, lower.tail = TRUE, log.p = FALSE){
  
  if (any(mu<=0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(theta<=0)) 
    stop(paste("theta must be positive", "\n", ""))
  if (any(nu<=0)) 
    stop(paste("nu must be positive", "\n", ""))
  
  if (log.p == TRUE) 
    p <- exp(p)
  else p <- p
  if (lower.tail == TRUE) 
    p <- p
  else p <- 1 - p
  if (any(p < 0) | any(p > 1)) 
    stop(paste("p must be between 0 and 1", "\n", ""))
  
  q <- {(1/mu)*log(1+(-1 + (1-p)^(-1))^(1/nu))}^(1/theta)
  q
}

#' @export
#' @rdname OW
rOW <- function(n,mu,theta,nu){
  if(any(n<=0))
    stop(paste("n must be positive","\n",""))
  if (any(mu<=0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(theta<=0)) 
    stop(paste("theta must be positive", "\n", ""))
  if (any(nu<=0)) 
    stop(paste("nu must be positive", "\n", ""))
  
  n <- ceiling(n)
  p <- runif(n)
  r <- qOW(p,mu,theta,nu)
  r
}
#' @export
#' @rdname OW
hOW<-function(x,mu,theta,nu){
  if (any(x<0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(mu<=0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(theta<=0)) 
    stop(paste("theta must be positive", "\n", ""))
  if (any(nu<=0)) 
    stop(paste("nu must be positive", "\n", ""))
  
  h <- dOW(x,mu,theta,nu, log = FALSE)/pOW(q=x,mu,theta,nu, lower.tail=FALSE, log.p = FALSE)
  h
}


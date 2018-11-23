#' @name SZMW
#' 
#' @title 
#' The Sarhan and Zaindins Modified Weibull Distribution
#' 
#' @description 
#' Density, distribution function, quantile function, 
#' random generation and hazard function for Sarhan and Zaindins modified weibull distribution with
#' parameters \code{mu}, \code{sigma} and \code{nu}.
#' 
#' @param x,q	vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. 
#' @param mu parameter one.    
#' @param sigma parameter two.
#' @param nu parameter three.
#' @param log,log.p	logical; if TRUE, probabilities p are given as log(p).	
#' @param lower.tail logical; if TRUE (default), probabilities are 
#' P[X <= x], otherwise, P[X > x].
#' 
#' @details 
#' The Sarhan and Zaindins modified weibull with parameters \code{mu}, 
#' \code{sigma} and \code{nu} has density given by
#' 
#' f(x)=(mu+sigma*nu*x^(nu-1))*exp(-mu*x-sigma*x^nu)
#' 
#' for x>0.
#' 
#' @return 
#' \code{dSZMW} gives the density, \code{pSZMW} gives the distribution 
#' function, \code{qSZMW} gives the quantile function, \code{rSZMW}
#' generates random deviates and \code{hSZMW} gives the hazard function.
#'  
#' @export
#' @examples 
#' 
#' ## The probability density function
#' curve(dSZMW(x, mu = 2, sigma = 1.5, nu = 0.2), from = 0, to = 2, ylim = c(0, 1.7), col = "red", las = 1, ylab = "The probability density function")
#' 
#' ## The cumulative distribution and the Reliability function
#' par(mfrow = c(1, 2))
#' curve(pSZMW(x, mu = 2, sigma = 1.5, nu = 0.2), from = 0, to = 2, ylim = c(0, 1), col = "red", las = 1, ylab = "The cumulative distribution function")
#' curve(pSZMW(x, mu = 2, sigma = 1.5, nu = 0.2, lower.tail = FALSE), from = 0, to = 2, ylim = c(0, 1), col = "red", las = 1, ylab = "The Reliability function")
#' 
#' ## The quantile function
#' p <- seq(from = 0, to = 0.99999, length.out = 100)
#' plot(x = qSZMW(p = p, mu = 2, sigma = 1.5, nu = 0.2), y = p, xlab = "Quantile", las = 1, ylab = "Probability")
#' curve(pSZMW(x, mu = 2, sigma = 1.5, nu = 0.2), from = 0, add = TRUE, col = "red")
#' 
#' ## The random function
#' hist(rSZMW(n = 1000, mu = 2, sigma = 1.5, nu = 0.2), freq = FALSE, xlab = "x", las = 1, main = "")
#' curve(dSZMW(x, mu = 2, sigma = 1.5, nu = 0.2),  from = 0, add = TRUE, col = "red")
#' 
#' ## The Hazard function
#' curve(hSZMW(x, mu = 2, sigma = 1.5, nu = 0.2), from = 0, to = 3, ylim = c(0, 8), col = "red", ylab = "The hazard function", las = 1)
#' 
SZMW <- function (mu.link = "log", sigma.link = "log", nu.link = "log") 
{
  mstats <- checklink("mu.link",    "Sarhan and Zaindins Modified Weibull", substitute(mu.link),    c("log", "own"))
  dstats <- checklink("sigma.link", "Sarhan and Zaindins Modified Weibull", substitute(sigma.link), c("log", "own"))
  vstats <- checklink("nu.link",    "Sarhan and Zaindins Modified Weibull", substitute(nu.link),    c("log", "own"))
  
  structure(list(family = c("SZMW", "Sarhan and Zaindins Modified Weibull"), 
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
                   exp1 <- mu + sigma*nu * y ^ (nu-1)
                   dexp1dm <- 1
                   dldm <- (1/exp1) * dexp1dm -y
                   dldm  
                 },
                 
                 d2ldm2 = function(y, mu, sigma, nu) {
                   exp1 <- mu + sigma*nu * y ^ (nu-1)
                   dexp1dm <- 1
                   d2ldm2 <- -(-((dexp1dm)^2/exp1 ^2))^2
                   d2ldm2
                 }, 
                 
                 dldd = function(y, mu, sigma, nu) {
                   exp1 <- mu + sigma*nu * y ^ (nu-1)
                   dexp1dd  <- nu*y^(nu-1)
                   dldd <-(1/exp1) * dexp1dd - y^nu
                   dldd
                 },
                 
                 d2ldd2 = function(y, mu, sigma, nu) {
                   exp1 <- mu + sigma*nu * y ^ (nu-1)
                   dexp1dd  <- nu*y^(nu-1) 
                   d2ldd2  <- -(-((dexp1dd)^2/exp1^2)-y^nu*log(y))^2
                   d2ldd2
                 }, 
                 
                 dldv = function(y, mu, sigma, nu) {
                   exp1 <- mu + sigma*nu * y ^ (nu-1) 
                   dexp1dv <- sigma * y^(nu-1)*(1+nu*log(y))
                   dldv <- (1/exp1)*dexp1dv-sigma*y^nu*log(y)
                   dldv
                 }, 
                 
                 d2ldv2 = function(y, mu, sigma, nu) {
                   exp1 <- mu + sigma*nu * y ^ (nu-1) 
                   dexp1dv <- sigma * y^(nu-1)*(1+nu*log(y)) 
                   d2exp1dv2  <- sigma*y^(nu-1)*log(y)*(nu*log(y)+2)
                   d2ldv2 <- ((exp1*d2exp1dv2-dexp1dv)/(exp1^2))-sigma*y^nu*(log(y))^2 
                   d2ldv2
                 }, 
                 
                 d2ldmdd = function(y, mu, sigma, nu) {
                   exp1 <- mu + sigma*nu * y ^ (nu-1)
                   dexp1dm <- 1
                   dexp1dd  <- nu*y^(nu-1)
                   d2ldmdd <- (dexp1dm*dexp1dd)/exp1^2
                   d2ldmdd
                 }, 
                 
                 d2ldmdv = function(y, mu, sigma, nu) {
                   exp1 <- mu + sigma*nu * y ^ (nu-1)
                   dexp1dm <- 1
                   dexp1dv <- sigma * y^(nu-1)*(1+nu*log(y))
                   d2ldmdv <- (dexp1dm*dexp1dv)/exp1^2
                   d2ldmdv
                 }, 
                 
                 d2ldddv = function(y, mu, sigma, nu) {
                   exp1 <- mu + sigma*nu * y ^ (nu-1)
                   dexp1dd  <- nu*y^(nu-1)
                   dexp1dv <- sigma * y^(nu-1)*(1+nu*log(y))
                   d2exp1dddv <- y^(nu-1)*(1+nu*log(y))
                   d2ldddv <- ((exp1*d2exp1dddv-dexp1dd*dexp1dv)/exp1^2)-y^nu*log(y)
                   d2ldddv
                 }, 
                 
                 G.dev.incr = function(y, mu, sigma, nu, ...) -2*dSZMW(y, mu, sigma, nu, log = TRUE), 
                 rqres = expression(rqres(pfun = "pSZMW", type = "Continuous",  y = y, mu = mu, sigma = sigma, nu = nu)), 
                 
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
#' @rdname SZMW
dSZMW<-function(x,mu,sigma,nu, log = FALSE){
  if (any(x<0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(mu<=0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma<=0)) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu<=0)) 
    stop(paste("nu must be positive", "\n", ""))
  
  loglik<- log(mu + sigma*nu*x^(nu-1)) - mu*x - sigma*x^nu
  
  if (log == FALSE) 
    density<- exp(loglik)
  else 
    density <- loglik
  return(density)
}
#' @export
#' @rdname SZMW
pSZMW <- function(q,mu,sigma,nu, lower.tail=TRUE, log.p = FALSE){
  if (any(q<0)) 
    stop(paste("q must be positive", "\n", ""))
  if (any(mu<=0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma<=0)) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu<=0)) 
    stop(paste("nu must be positive", "\n", ""))
  
  cdf <- 1- exp(-mu*q -sigma*(q^nu))
  if (lower.tail == TRUE) 
    cdf <- cdf
  else cdf <- 1 - cdf
  if (log.p == FALSE) 
    cdf <- cdf
  else cdf <- log(cdf)
  cdf
  
}
#' @export
#' @rdname SZMW
qSZMW <- function(p, mu,sigma,nu, lower.tail = TRUE, log.p = FALSE) {
  if (any(mu<=0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma<=0)) 
    stop(paste("sigma must be positive", "\n", ""))
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
  
  fda <- function(x,mu,sigma,nu){
    1- exp(-mu*x - sigma*(x^nu))
  }
  
  fda1 <- function(x, mu,sigma,nu, p) {fda(x, mu,sigma,nu) - p}
  
  r_de_la_funcion <- function(mu,sigma,nu, p) {
    uniroot(fda1, interval=c(0,1e+06), mu,sigma,nu, p)$root
  }
  
  r_de_la_funcion <- Vectorize(r_de_la_funcion)
  q <- r_de_la_funcion(mu,sigma,nu, p)
  q
  
}
#' @export
#' @rdname SZMW
rSZMW<- function(n,mu,sigma,nu){
  if (any(mu<=0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma<=0)) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu<=0)) 
    stop(paste("nu must be positive", "\n", ""))
  
  n <- ceiling(n)
  p <- runif(n)
  r <- qSZMW(p, mu,sigma,nu)
  r
}
#' @export
#' @rdname SZMW
hSZMW<-function(x,mu,sigma,nu){
  if (any(x<0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(mu<=0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma<=0)) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu<=0)) 
    stop(paste("nu must be positive", "\n", ""))
  
  h <- dSZMW(x,mu,sigma,nu, log = FALSE)/pSZMW(q=x,mu,sigma,nu, lower.tail=FALSE, log.p = FALSE)
  h  
}



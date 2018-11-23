#' @name MWEx
#' 
#' @title 
#' The Modified Weibull Extension Distribution 
#' 
#' @description 
#' Density, distribution function, quantile function, 
#' random generation and hazard function for the modified weibull extension distribution with
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
#' The modified weibull extension distribution with parameters \code{mu}, \code{sigma}
#' and \code{nu} has density given by
#' 
#' f(x) = nu*sigma*((x/mu)^(sigma-1))*exp((x/mu)^sigma + nu*mu*(1-exp((x/mu)^sigma)))
#' 
#' for x > 0.
#' 
#' @return 
#' \code{dMWEx} gives the density, \code{pMWEx} gives the distribution 
#' function, \code{qMWEx} gives the quantile function, \code{rMWEx}
#' generates random deviates and \code{hMWEx} gives the hazard function.
#' 
#' @export
#' @examples  
#' ## The probability density function  
#' curve(dMWEx(x, mu = 1/0.5, sigma = 3, nu = 2),
#'       from = 0, to = 2.5, ylim = c(0, 1.5), 
#'       col = "red", las = 1, 
#'       ylab = "The probability density function")
#' 
#' ## The cumulative distribution and the Reliability function
#' par(mfrow = c(1, 2))
#' curve(pMWEx(x, mu = 1/0.5, sigma = 3, nu = 2), 
#'       from = 0, to = 2.5, ylim = c(0, 1), 
#'       col = "red", las = 1, 
#'       ylab = "The cumulative distribution function")
#' curve(pMWEx(x, mu = 1/0.5, sigma = 3, nu = 2,  lower.tail = FALSE), 
#'       from = 0, to = 2.5, ylim = c(0, 1), 
#'       col = "red", las = 1, ylab = "The Reliability function")
#' 
#' ## The quantile function
#' p <- seq(from = 0, to = 0.998, length.out = 100)
#' plot(x = qMWEx(p, mu = 1/0.5, sigma = 3, nu = 2), 
#'      y = p, xlab = "Quantile", las = 1, 
#'      ylab = "Probability")
#' curve(pMWEx(x, mu = 1/0.5, sigma = 3, nu = 2), 
#'       from = 0, add = TRUE, col = "red")
#' 
#' ## The random function
#' hist(rMWEx(n = 10000, mu = 1/0.5, sigma = 3, nu = 2), 
#'      freq = FALSE, ylim = c(0, 1.5), 
#'      xlab = "x", las = 1, main = "")
#' curve(dMWEx(x, mu = 1/0.5, sigma = 3, nu = 2),  
#'       from = 0, ylim = c(0, 2.5), add = T, col = "red")
#' 
#' ## The Hazard function
#' curve(hMWEx(x, mu = 1/0.5, sigma = 3, nu = 2), from = 0, to = 1.7, ylim = c(0, 12), col = "red", ylab = "The hazard function", las = 1)
MWEx <- function (mu.link = "log", sigma.link = "log", nu.link = "log") 
{
  mstats <- checklink("mu.link",    "Modified Weibull Extension", substitute(mu.link),    c("log", "own"))
  dstats <- checklink("sigma.link", "Modified Weibull Extension", substitute(sigma.link), c("log", "own"))
  vstats <- checklink("nu.link",    "Modified Weibull Extension", substitute(nu.link),    c("log", "own"))
  
  structure(list(family = c("MWEx", "Modified Weibull Extension"), 
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
                   exp1 <- (y/sigma)^nu
                   exp2 <- mu*sigma*(1-exp(exp1)) 
                   dexp2dm <- sigma*(1-exp(exp1))
                   dldm <- (1/mu) + dexp2dm
                   dldm
                 },
                 
                 d2ldm2 = function(y, mu, sigma, nu) {
                   d2ldm2 <-(1/mu^2)   
                   d2ldm2  
                 }, 
                 
                 dldd = function(y, mu, sigma, nu) {
                   exp1 <- (y/sigma)^nu
                   exp2 <- mu*sigma*(1-exp(exp1))
                   dexp1dd <- -((v^y)/sigma^2)*(y/sigma)^(v-1)
                   dexp2dd <-  mu*((1-exp(exp1))-sigma*exp(exp1)*dexp1dd)
                   dldd <- (nu-1/sigma)+dexp1dd+dexp2dd
                   dldd
                 },
                 
                 d2ldd2 = function(y, mu, sigma, nu) {
                   exp1 <- (y/sigma)^nu
                   exp2 <- mu*sigma*(1-exp(exp1))
                   dexp1dd <- -((v^y)/sigma^2)*(y/sigma)^(v-1)
                   dexp2dd <-  mu*((1-exp(exp1))-sigma*exp(exp1)*dexp1dd)
                   d2exp1dd2 <- nu*y*((2/sigma^3)*(y/sigma)^(nu-1)+((y*(nu-1))/sigma^4)*(y/sigma)^(nu-2))
                   d2exp2dd2 <- -(-mu*exp(exp1)*((dexp1dd)+sigma*(dexp1dd)^2+sigma*d2exp1dd2))^2
                   d2ldd2 <- ((nu-1)/sigma^2)+d2exp1dd2+ d2exp2dd2
                   d2ldd2
                 }, 
                 
                 dldv = function(y, mu, sigma, nu) {
                   exp1 <- (y/sigma)^nu
                   exp2 <- mu*sigma*(1-exp(exp1))
                   dexp1dv <- (y/sigma)^nu*log10(y/sigma)
                   dexp2dv<- -mu*sigma*exp(exp1)* dexp1dv  
                   dldv <- 1/nu+ log(y/sigma)+dexp1dv+dexp2dv
                   dldv
                 }, 
                 
                 d2ldv2 = function(y, mu, sigma, nu) {
                   exp1 <- (y/sigma)^nu
                   exp2 <- mu*sigma*(1-exp(exp1))
                   dexp1dv <- (y/sigma)^nu*log10(y/sigma)
                   d2exp1dv2 <- (y/sigma)^nu*(log10(y/sigma))^2
                   d2exp2dv2 <- -mu*sigma*exp(exp1)*((dexp1dv)^2+d2exp1dv2)
                   d2ldv2 <- -(-1/nu^2+d2exp1dv2+d2exp2dv2)^2
                   d2ldv2
                 }, 
                 
                 d2ldmdd = function(y, mu, sigma, nu) {
                   exp1 <- (y/sigma)^nu
                   exp2 <- mu*sigma*(1-exp(exp1))
                   dexp1dd <- -((v^y)/sigma^2)*(y/sigma)^(v-1)   
                   d2exp2dmdd  <- 1-exp(exp1)*(1+sigma *dexp1dd)
                   d2ldmdd  <-d2exp2dmdd
                   d2ldmdd
                 }, 
                 
                 d2ldmdv = function(y, mu, sigma, nu) {
                   exp1 <- (y/sigma)^nu
                   exp2 <- mu*sigma*(1-exp(exp1)) 
                   dexp1dv <- (y/sigma)^nu*log(y/sigma)
                   d2exp2dmdv  <- -sigma*exp(exp1) * dexp1dv
                   d2ldmdv <- d2exp2dmdv
                   d2ldmdv
                 }, 
                 
                 d2ldddv = function(y, mu, sigma, nu) {
                   exp1 <- (y/sigma)^nu
                   exp2 <- mu*sigma*(1-exp(exp1)) 
                   dexp1dv <- (y/sigma)^nu*log10(y/sigma)
                   dexp1dd <- -((v^y)/sigma^2)*(y/sigma)^(v-1)
                   d2exp1dddv <- -(y/sigma^2)*(y/sigma)^(nu-1)*(1+nu*log(y/sigma))
                   d2exp2dddv <- -mu*exp(exp1)*dexp1dv-sigma*exp(exp1)*(dexp1dv*dexp1dd*d2exp1dddv)
                   d2ldddv <- -(-1/sigma+d2exp1dddv+d2exp2dddv)^2
                   d2ldddv
                 }, 
                 
                 G.dev.incr = function(y, mu, sigma, nu, ...) -2*dMWEx(y, mu, sigma, nu, log = TRUE), 
                 rqres = expression(rqres(pfun = "pMWEx", type = "Continuous",  y = y, mu = mu, sigma = sigma, nu = nu)), 
                 
                 mu.initial = expression( mu <-  rep(0.5, length(y)) ), 
                 sigma.initial = expression( sigma <- rep(0.5, length(y)) ), 
                 nu.initial = expression( nu <- rep(0.5, length(y)) ), 
                 
                 mu.valid = function(mu) all(mu >  0), 
                 sigma.valid = function(sigma) all(sigma >  0), 
                 nu.valid = function(nu) all(nu > 0), 
                 
                 y.valid = function(y) all(y >= 0)
  ), 
  class = c("gamlss.family", "family"))
}
#' @export
#' @rdname MWEx
dMWEx<-function(x,mu,sigma,nu, log = FALSE){
  if (any(x<0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(mu<=0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma<=0)) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu<=0)) 
    stop(paste("nu must be positive", "\n", ""))
  
  
  loglik<-log(nu*sigma) + (sigma-1)*log(x/mu) +
    (x/mu)^sigma + nu*mu*(1-exp((x/mu)^sigma))
  
  if (log == FALSE) 
    density<- exp(loglik)
  else 
    density <- loglik
  return(density)
}

#' @export
#' @rdname MWEx

pMWEx<- function(q,mu,sigma,nu, lower.tail=TRUE, log.p = FALSE){
  if (any(q<0)) 
    stop(paste("q must be positive", "\n", ""))
  if (any(mu<=0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma<=0)) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu<=0)) 
    stop(paste("nu must be positive", "\n", ""))
  
  cdf <- 1- exp(nu*mu*(1-exp((q/mu)^sigma)))
  
  if (lower.tail == TRUE) 
    cdf <- cdf
  else cdf <- 1 - cdf
  if (log.p == FALSE) 
    cdf <- cdf
  else cdf <- log(cdf)
  cdf
}

#' @export
#' @rdname MWEx
qMWEx <- function(p,mu,sigma,nu, lower.tail = TRUE, log.p = FALSE){
  
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
  
  q <- mu*(log(1-((1/(nu*mu))*log(1-p))))^(1/sigma)
  q
}


#' @export
#' @rdname MWEx
rMWEx <- function(n,mu,sigma,nu){
  if(any(n<=0))
    stop(paste("n must be positive","\n",""))
  if (any(mu<=0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma<=0)) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu<=0)) 
    stop(paste("nu must be positive", "\n", ""))
  n <- ceiling(n)
  p <- runif(n)
  r <- qMWEx(p, mu,sigma,nu)
  r
}

#' @export
#' @rdname MWEx
hMWEx<-function(x,mu,sigma,nu){
  if (any(x<0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(mu<=0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma<=0)) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu<=0)) 
    stop(paste("nu must be positive", "\n", ""))
  
  h <- dMWEx(x,mu,sigma,nu, log = FALSE)/pMWEx(q=x,mu,sigma,nu, lower.tail=FALSE, log.p = FALSE)
  h
}

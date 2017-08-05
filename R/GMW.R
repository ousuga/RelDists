#' @name GMW
#' 
#' @title 
#' The Generalized modified Weibull Distribution 
#' 
#' @description 
#' Density, distribution function, quantile function, random generation and hazard function for the generalized 
#' modified weibull distribution with parameters \code{mu}, \code{sigma}, \code{nu} and \code{tau}.
#' 
#' @param x,q	vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. 
#' @param mu parameter one.
#' @param sigma parameter two.
#' @param nu parameter three.
#' @param tau parameter four.
#' @param log,log.p	logical; if TRUE, probabilities p are given as log(p).	
#' @param lower.tail logical; if TRUE (default), probabilities are 
#' P[X <= x], otherwise, P[X > x].
#' 
#' @details 
#' The generalized modified weibull with parameters \code{mu}, \code{sigma}, 
#' \code{nu} and \code{tau} has density given by
#' 
#' f(x)= mu*sigma*x^(nu-1)*(nu+tau*x)*exp(tau*x-mu*x^(nu)*exp(tau*x))*
#' (1-exp(mu*x^(nu)*exp(tau*x)))^(sigma-1)
#' 
#' for x>0.
#' 
#' #' @return 
#' \code{dGMW} gives the density, \code{pGMW} gives the distribution 
#' function, \code{qGMW} gives the quantile function, \code{rGMW}
#' generates random deviates and \code{hGMW} gives the hazard function.
#' 
#' @export
#' @examples  
#'## The probability density function
#' curve(dGMW(x, mu = 2, sigma = 0.5, nu = 2, tau = 1.5), from = 0, to = 0.8, ylim = c(0, 6), col = "red", las = 1, ylab = "The probability density function") 
#' 
#' ## The cumulative distribution and the Reliability function
#' par(mfrow = c(1, 2))
#' curve(pGMW(x, mu = 2, sigma = 0.5, nu = 2, tau = 1.5), from = 0, to = 1.2, col = "red", las = 1, ylab = "The cumulative distribution function")
#' curve(pGMW(x, mu = 2, sigma = 0.5, nu = 2, tau = 1.5, lower.tail = FALSE), from = 0, to = 1.2, col = "red", las = 1, ylab = "The Reliability function")

#' ## The quantile function
#' p <- seq(from = 0, to = 0.99999, length.out = 100)
#' plot(x = qGMW(p, mu = 2, sigma = 0.5, nu = 2, tau = 0.3), y = p, xlab = "Quantile", las = 1, ylab = "Probability")
#' curve(pGMW(x, mu = 2, sigma = 0.5, nu = 2, tau = 0.3),  from = 0, add = TRUE, col="red")
#' 
#' ## The random function
#' hist(rGMW(n = 1000, mu = 2, sigma = 0.5, nu = 2,tau = 0.3), freq = FALSE, xlab = "x", main ="", las = 1)
#' curve(dGMW(x, mu = 2, sigma = 0.5, nu = 2, tau = 0.3),  from = 0, add = TRUE, col = "red")
#' 
#' ## The Hazard function
#'curve(hGMW(x, mu = 2, sigma = 1.5, nu = 2, tau = 0.8), from = 0, to = 1, ylim = c(0, 16), col = "red", ylab = "The Hazard function", las = 1)
#'
GMW <- function (mu.link = "log", sigma.link = "log", nu.link = "log", tau.link = "log") 
{
  mstats <- checklink("mu.link",    "Generalized Modified Weibull", substitute(mu.link),    c("log", "own"))
  dstats <- checklink("sigma.link", "Generalized Modified Weibull", substitute(sigma.link), c("log", "own"))
  vstats <- checklink("nu.link",    "Generalized Modified Weibull", substitute(nu.link),    c("log", "own"))
  tstats <- checklink("tau.link",   "Generalized Modified Weibull", substitute(tau.link),    c("log", "own"))
  
  structure(list(family = c("GMW", "Generalized Modified Weibull"), 
                 parameters = list(mu = TRUE, sigma = TRUE, nu = TRUE, tau = TRUE), 
                 nopar = 4, 
                 type = "Continuous", 
                 
                 mu.link    = as.character(substitute(mu.link)), 
                 sigma.link = as.character(substitute(sigma.link)), 
                 nu.link    = as.character(substitute(nu.link)), 
                 tau.link   = as.character(substitute(tau.link)), 
                 
                 mu.linkfun    = mstats$linkfun, 
                 sigma.linkfun = dstats$linkfun, 
                 nu.linkfun    = vstats$linkfun, 
                 tau.linkfun   = tstats$linkfun, 
                 
                 mu.linkinv    = mstats$linkinv, 
                 sigma.linkinv = dstats$linkinv, 
                 nu.linkinv    = vstats$linkinv, 
                 tau.linkinv   = tstats$linkinv,
                 
                 mu.dr    = mstats$mu.eta, 
                 sigma.dr = dstats$mu.eta, 
                 nu.dr    = vstats$mu.eta, 
                 tau.dr   = tstats$mu.eta,
                 
                 
                 dldm = function(y, mu, sigma, nu, tau) {
                   dldm<- (((sigma-1)*exp(tau+y)*(y^nu))/(exp(mu*exp(tau+y)*(y^nu))-1))+(1/(mu-sigma))-exp(tau)*(y^nu)
                   dldm
                 },
                 
                 d2ldm2 = function(y, mu, sigma, nu, tau) {
                   d2ldm2<- -(((sigma-1)*y^(2*nu)*exp(mu*exp(tau+y)*(y^nu)+2*tau+2*y))/(exp(x*exp(tau+y)*(y^nu))-1)^2)-(1/(sigma-mu)^2)
                   d2ldm2
                 }, 
                 
                 dldd = function(y, mu, sigma, nu, tau) {
                   dldd<- (1/(sigma-mu))+log(1-exp(x*(-exp(tau+y)*(y^sigma))))
                   dldd
                 },
                 
                 d2ldd2 = function(y, mu, sigma, nu, tau) {
                   d2ldd2<- -(1/((sigma-mu)^2))
                   d2ldd2
                 }, 
                 
                 dldv = function(y, mu, sigma, nu, tau) {
                   dldv<- log(y)*((((sigma-1)*mu*exp(tau+y)*(y^nu))/(exp(mu*exp(tau+y)*(y^nu))-1)) -exp(tau)*mu*(y^nu)+1)+(1/tau*y+nu) 
                   dldv
                 }, 
                 
                 d2ldv2 = function(y, mu, sigma, nu, tau) {
                   d2ldv2<- -(((sigma-1)*mu^2*y(2*nu)*log^2(y)*exp(mu*exp(tau+y)*(y^nu)+2*tau+2*y))/(exp(mu*exp(tau+y)*(y^nu))-1))+((exp(tau)*mu*(y^nu)*log^2(y)*(sigma*exp(y))-exp(mu*exp(tau+y)*(y^nu))-exp(y)+1)/(exp(mu*exp(tau+y)*(y^nu))-1))-(1/(tau*y+nu)^2)
                   d2ldv2
                 }, 
                 
                 dldt = function(y, mu, sigma, nu, tau) {
                   dldt<- (((nu-1)*mu*exp(tau+y)*(y^nu))/(exp(mu*exp(tau+y)*y^nu)-1))-exp(tau)*mu*(y^nu)+(y/(tau*y+nu))+y
                   dldt
                 },
                 
                 d2ldt2 = function(y, mu, sigma, nu, tau) {
                   d2ldt2<- -(((nu-1)*mu*exp(tau+y)*(y^nu)*(mu*(y^nu)*exp(mu*exp(tau+y)*(y^nu)+tau+y))-exp(x*exp(tau+y)*(y^nu))+1)/((exp(x*exp(tau6y)*(y^nu))-1)^2)) - exp(tau)*mu*(y^nu)- ((y^2)/(tau*y+nu)^2)   
                 },
                 
                 d2ldmdd = function(y, mu, sigma, nu, tau) {
                   d2ldmdd<- 1/(mu-sigma)^2+((exp(tau+y))*(y^nu))/((exp(mu*exp(tau+y)))*(y^(nu))-1)
                   d2ldmdd     
                 }, 
                 
                 d2ldmdv = function(y, mu, sigma, nu, tau) {
                   d2ldmdv<- (((sigma-1)*(exp(tau+y))*(y^nu)*log(y))/(exp(mu*exp(tau+y)))*(y^nu)-1)-(((sigma-1)*mu*(y^(2*nu))*log(y)*exp(mu*(exp(tau+y))*(y^nu)+2*tau+2*y))/((exp(mu*(exp(tau+y))*(y^nu))-1)^2))-(exp(tau)*(y^nu)*log(y)) 
                   d2ldmdv  
                 }, 
                 
                 d2ldmdt = function(y, mu, sigma, nu, tau) {
                   d2ldmdt<- (((sigma-1)*(exp(tau+y))*(y^nu))/(exp(mu*(exp(tau+y))*(y^nu))-1))-(((sigma-1)*mu*(y^2*nu)*exp(mu*(exp(tau+y))*(y^nu)+2*tau+2*y))/(exp(mu*(exp(tau+y))*(y^nu))-1))-(exp(tau)*(y^nu))
                   d2ldmdt
                 }, 
                 
                 d2ldddv = function(y, mu, sigma, nu, tau) {
                   d2ldddv<- 0
                   d2ldddv
                 }, 
                 d2ldddt = function(y, mu, sigma, nu, tau) {
                   d2ldddt<- (mu*(y^sigma)*exp(tau+y))/((exp(mu*(y^sigma))*exp(tau+y))-1)
                   d2ldddt
                 }, 
                 
                 d2ldvdt = function(y, mu, sigma, nu, tau) {
                   d2ldvdt<- log(y)*(-(((sigma-1)*(mu^2)*(y^2*nu)*(exp(mu*exp(tau+y))*(y^nu)*2*tau+2*y))/(((exp(mu*exp(tau+y))*(y^nu))-1)^2))+(((sigma-1)*mu*(exp(tau+y))*(y^nu))/((exp(mu*(exp(tau+y))*(y^nu))-1)))-(exp(tau))*mu*(y^nu))-(y/(mu^2))
                   d2ldvdt
                 }, 
                 
                 G.dev.incr = function(y, mu, sigma, nu, tau, ...) -2*dGMW(y, mu, sigma, nu, tau, log = TRUE), 
                 rqres = expression(rqres(pfun = "pGMW", type = "Continuous",  y = y, mu = mu, sigma = sigma, nu = nu, tau = tau)), 
                 
                 mu.initial = expression( mu <-  rep(0.5, length(y)) ), 
                 sigma.initial = expression( sigma <- rep(0.5, length(y)) ), 
                 nu.initial = expression( nu <- rep(0.5, length(y)) ), 
                 tau.initial = expression( tau <- rep(0.5, length(y)) ),
                 
                 mu.valid = function(mu) all(mu >  0), 
                 sigma.valid = function(sigma) all(sigma >  0), 
                 nu.valid = function(nu) all(nu >= 0), 
                 tau.valid = function(tau) all(tau >= 0), 
                 
                 y.valid = function(y) all(y > 0)
  ), 
  class = c("gamlss.family", "family"))
}

#' @export
#' @rdname GMW
dGMW<-function(y,mu,sigma,nu,tau, log = FALSE){
  if (any(y<0)) 
    stop(paste("y must be positive", "\n", ""))
  if (any(mu<=0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma<=0)) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu<0)) 
    stop(paste("nu must be positive", "\n", ""))
  if (any(tau<0)) 
    stop(paste("tau must be positive", "\n", ""))
  
  log_fy<-log(mu*sigma) + (nu-1)*log(y) + log(nu +tau*y) +
    tau*y - mu*(y^nu)*exp(tau*y)  +
    (sigma-1)*log(1-exp(-mu*(y^nu)*exp(tau*y) ))
  
  if (log == FALSE) 
    density<- exp(log_fy)
  else 
    density <- log_fy
  return(density)
}

#' @export
#' @rdname GMW
pGMW <- function(q,mu,sigma,nu,tau, lower.tail=TRUE, log.p = FALSE){
  if (any(q<0)) 
    stop(paste("q must be positive", "\n", ""))
  if (any(mu<=0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma<=0)) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu<0)) 
    stop(paste("nu must be positive", "\n", ""))
  if (any(tau<0)) 
    stop(paste("tau must be positive", "\n", ""))
  
  cdf <- (1-exp(-mu*(q^nu)*exp(tau*q)))^sigma
  if (lower.tail == TRUE) 
    cdf <- cdf
  else cdf <- 1 - cdf
  if (log.p == FALSE) 
    cdf <- cdf
  else cdf <- log(cdf)
  cdf
}


#' @export
#' @rdname GMW
qGMW <- function(p, mu,sigma,nu,tau, lower.tail = TRUE, log.p = FALSE) {
  if (any(mu<=0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma<=0)) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu<0)) 
    stop(paste("nu must be positive", "\n", ""))
  if (any(tau<0)) 
    stop(paste("tau must be positive", "\n", ""))
  
  if (log.p == TRUE) 
    p <- exp(p)
  else p <- p
  if (lower.tail == TRUE) 
    p <- p
  else p <- 1 - p
  if (any(p < 0) | any(p > 1)) 
    stop(paste("p must be between 0 and 1", "\n", ""))
  
  fda <- function(x,mu,sigma,nu,tau){
    (1-exp((-mu*(x^nu))*exp(tau*x)))^sigma
  }
  
  fda1 <- function(x, mu,sigma,nu,tau, p) {fda(x, mu,sigma,nu,tau) - p}
  
  r_de_la_funcion <- function(mu,sigma,nu,tau, p) {
    uniroot(fda1, interval=c(0,1e+06), mu,sigma,nu,tau, p)$root
  }
  
  r_de_la_funcion <- Vectorize(r_de_la_funcion)
  q <- r_de_la_funcion(mu,sigma,nu,tau, p)
  q
  
}


#' @export
#' @rdname GMW
rGMW <- function(n,mu,sigma,nu,tau){
  if (any(mu<=0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma<=0)) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu<0)) 
    stop(paste("nu must be positive", "\n", ""))
  if (any(tau<0)) 
    stop(paste("tau must be positive", "\n", ""))
  
  n <- ceiling(n)
  p <- runif(n)
  r <- qGMW(p, mu,sigma,nu,tau)
  r
}


#' @export
#' @rdname GMW
hGMW<-function(x,mu,sigma,nu,tau, log = FALSE){
  if (any(x<0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(mu<=0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma<=0)) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu<0)) 
    stop(paste("nu must be positive", "\n", ""))
  if (any(tau<0)) 
    stop(paste("tau must be positive", "\n", ""))
  
  h <- dGMW(x,mu,sigma,nu,tau, log = FALSE)/pGMW(q=x,mu,sigma,nu,tau, lower.tail=FALSE, log.p = FALSE)
  h
}




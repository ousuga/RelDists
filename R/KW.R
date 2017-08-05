#' @name KW
#' 
#' @title 
#' The Kumaraswamy Weibull Distribution 
#' 
#' @description 
#' Density, distribution function, quantile function, 
#' random generation and hazard function for the kumaraswamy weibull distribution with
#' parameters \code{mu}, \code{sigma}, \code{nu} and \code{tau}.
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
#' The kumaraswamy weibull with parameters \code{mu}, \code{sigma}, \code{nu} and \code{tau}
#' has density given by
#' 
#' f(x)=nu*tau*mu*sigma*x^(sigma-1)*exp(-mu*x^sigma)*(1-exp(-mu*x^sigma))^(nu-1)*
#' (1-(1-exp(-aplha*x^sigma))^nu)^(tau-1)
#' 
#' for x>0.
#' 
#' @return 
#' \code{dKW} gives the density, \code{pKW} gives the distribution 
#' function, \code{qKW} gives the quantile function, \code{rKW}
#' generates random deviates and \code{hKW} gives the hazard function.
#' 
#' @export
#' @examples 
#' ## The probability density function  
#' curve(dKW(x, mu = 3, sigma = 0.8, nu = 2, tau = 1.5), from = 0, to = 3, ylim = c(0, 2.5), col = "red", las = 1, ylab = "The probability density function")
#'
#' ## The cumulative distribution and the Reliability function
#' par(mfrow = c(1, 2))
#' curve(pKW(x, mu = 3, sigma = 0.8, nu = 2, tau = 1.5), from = 0, to = 3, ylim = c(0, 1), col = "red", las = 1, ylab = "The cumulative distribution function")
#' curve(pKW(x, mu = 3, sigma = 0.8, nu = 2, tau = 1.5, lower.tail = FALSE), from = 0, to = 3, ylim = c(0, 1), col = "red", las = 1, ylab = "The Reliability function")
#' 
#' ## The quantile function
#' p <- seq(0,0.99999, length.out=100)
#' plot(x=qKW(p, mu = 3, sigma = 0.8, nu = 2, tau = 1.5), y = p, xlab = "Quantile", las = 1, ylab = "Probability")
#' curve(pKW(x, mu = 3, sigma = 0.8, nu = 2, tau = 1.5), from = 0, add = TRUE, col = "red")
#' 
#' ## The random function
#' hist(rKW(10000, mu = 3, sigma = 0.8, nu = 2, tau = 1.5), freq = FALSE, xlab = "x", las = 1, ylim = c(0, 2.5), main = "")
#' curve(dKW(x, mu = 3, sigma = 0.8, nu = 2, tau = 1.5),  from = 0, add = TRUE, col = "red" )
#' 
#' ## The Hazard function
#' curve(hKW(x, mu = 3, sigma = 0.8, nu = 2, tau = 1.5), from = 0, to = 2.1, ylim = c(0, 4), col = "red", ylab = "The Hazard function", las = 1)
KW <- function (mu.link = "log", sigma.link = "log", nu.link = "log", tau.link = "log") 
{
  mstats <- checklink("mu.link",    "Kumaraswamy Weibull", substitute(mu.link),    c("log", "own"))
  dstats <- checklink("sigma.link", "Kumaraswamy Weibull", substitute(sigma.link), c("log", "own"))
  vstats <- checklink("nu.link",    "Kumaraswamy Weibull", substitute(nu.link),    c("log", "own"))
  tstats <- checklink("tau.link",   "Kumaraswamy Weibull", substitute(tau.link),   c("log", "own"))
  
  structure(list(family = c("KW", "Kumaraswamy Weibull"), 
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
                   dldm<- -(((sigma-1)*(1-exp(-nu*(y^tau)))^mu *log(1-exp(-nu*(y^tau))))/(1-(1-exp(-nu*(y^tau)))))+log(1-exp(-nu*(y^tau)))+(1/mu)
                   dldm
                 },
                 
                 d2ldm2 = function(y, mu, sigma, nu, tau) {
                   d2ldm2<- -(((sigma-1)*(1-exp(-nu*(y^tau))*log^2(1-exp(-nu*(y^tau)))))/((1-exp(-nu*(y^tau)))^mu)^2)-(1/(mu^2))
                   d2ldm2
                 }, 
                 
                 dldd = function(y, mu, sigma, nu, tau) {
                   dldd<- (sigma*log(1-(1-exp(-nu*(y^tau)))^mu)+1)/sigma  
                   dldd
                 },
                 
                 d2ldd2 = function(y, mu, sigma, nu, tau) {
                   d2ldd2<- -(1/sigma^2)
                   d2ldd2
                 }, 
                 
                 dldv = function(y, mu, sigma, nu, tau) {
                   dldv<- -(((sigma-1)*mu*(y^tau)*exp(-nu*(y^tau))*(1-exp(-nu*(y^tau)))^(mu-1))/(1-(1-exp(-nu*(y^tau)))^mu))+(((mu-1)*(y^tau)*exp(-nu*(y^tau)))/(1-exp(-nu*(y^tau))))-(y^tau)+(1/nu)
                   dldv
                 }, 
                 
                 d2ldv2 = function(y, mu, sigma, nu, tau) {
                   d2ldv2<- (sigma-1)*(((-mu^2)*y^(2*tau)*exp(-2*nu*(y^tau))*(1-exp(-nu*(y^tau)))^(2*mu-2))/((1-(1-exp(-nu*(y^tau)))^mu)^2)*((mu-1)*mu*(y^(2*tau))*exp(-2*nu*(y^tau))*(1-exp(-nu*(y^tau)))^(mu-2))/(1-(1-exp(-nu*(y^tau))))+((mu*y^(2*tau))*exp(-nu*(y^tau))*(1-exp(-nu*(y^tau)))^(mu-1))/(1-(1-exp(-nu*(y^tau)))^mu))+(mu-1)*(-((y^(2*tau))*exp(-nu*(y^tau)))/(1-exp(-nu*(y^tau)))-(((-y^(2*tau))*exp(-2*nu*(y^tau)))/((1-exp(-nu*(y^tau))))))-(1/(nu^2))
                   d2ldv2
                 }, 
                 
                 dldt = function(y, mu, sigma, nu, tau) {
                   dldt<- -(((sigma-1)*nu*mu*(y^tau)*log(y)*exp(-nu*(y^tau))*(1-exp(-nu*(y^tau)))^(mu-1))/(1-(1-exp(-nu*(y^tau)))^mu))
                   dldt
                 },
                 
                 d2ldt2 = function(y, mu, sigma, nu, tau) {
                   d2ldt2<- (1/(tau^2))-((sigma-1)*nu*mu*(y^tau)*log^2(y)*((1-exp(-nu*(y^tau)))^mu)*(nu*mu)*(y^tau)+((1-exp(-nu*(y^tau)))^mu)+exp(nu*(y^tau))*(nu*(y^tau)-1)*(((1-exp(-nu*(y^tau)))^mu)-1)-1)/(((exp(nu*(y^tau))-1)^2)*(((1-exp(-nu*(y^tau)))^mu)-1)^2)-(((nu*(mu-1))*(y^tau)*log^2(y)*(exp(nu*(y^tau))*(nu*(y^tau))-1))/((exp(nu*(y^tau))-1)^2))-nu*(y^tau)*log^2(y)
                   d2ldt2
                 },
                 
                 d2ldmdd = function(y, mu, sigma, nu, tau) {
                   d2ldmdd<- -exp(nu*(y^tau))*(1-exp(-nu*(y^tau)))^mu*log(1-exp(-nu*(y^tau))) 
                   d2ldmdd     
                 }, 
                 
                 d2ldmdv = function(y, mu, sigma, nu, tau) {
                   d2ldmdv<- (sigma-1)*(y^tau)*(-(1-exp(-nu*(y^tau)))^(mu-1))-((sigma-1)*(y^tau)*exp(nu*(y^tau)))^mu*log(1-exp(-nu*(y^tau)))-(sigma-1)*mu*(y^tau)*(1-exp(-nu*(y^tau)))^(mu-1)*log(1-exp(-nu*(y^tau)))+ (((y^tau)*exp(-nu*(y^tau)))/(1-exp(-nu*(y^tau)))) 
                   d2ldmdv  
                 }, 
                 
                 d2ldmdt = function(y, mu, sigma, nu, tau) {
                   d2ldmdt<- (((sigma-1)*(exp(nu*(y^tau)))*(1-exp(-nu*(y^tau)))^mu)+(sigma-1)*((exp(nu*(y^tau)))+mu-1)*((1-exp(-nu*(y^tau)))^mu)*(log(1-exp(-nu*(y^tau)))-1))
                   d2ldmdt
                 }, 
                 
                 d2ldddv = function(y, mu, sigma, nu, tau) {
                   d2ldddv<- ((-mu*(y^tau)*exp(-nu*(y^tau))*(1-exp(-nu*(y^tau)))^mu-1)/(1-(1-exp(-nu*(y^tau)))^mu))
                   d2ldddv
                 }, 
                 
                 d2ldddt = function(y, mu, sigma, nu, tau) {
                   d2ldddt<- -((nu*mu*(y^tau)*log(y)*exp(-nu*(y^tau))*(1-exp(nu*(y^tau)))^mu-1)/(1-(1-exp(-nu*(y^tau)))^mu))
                   d2ldddt
                 }, 
                 
                 d2ldvdt = function(y, mu, sigma, nu, tau) {
                   d2ldvdt<- (-(y^tau)*log(y)*(exp(nu(y^tau)))*((1-exp(-nu*(y^tau))^mu)-1)*(nu*(y^tau)*(mu*(sigma*((1-exp(-nu*(y^tau)))^mu)-1)-((1-exp(-nu*(y^tau)))^mu)+1)-sigma*mu*(1-exp(-nu*(y^tau)))^mu-((1-exp(-nu*(y^tau)))^mu)+mu+1)+mu*(sigma*(1-exp(-nu*(y^tau)))^mu(nu*mu*(y^tau)+(1-exp(-nu*(y^tau)))^mu+1)+exp(2*nu*(y^tau))*((1-exp(-nu*(y^tau)))^mu)-1)^2)/(((exp(nu*(y^tau))-1)^2)*((1-exp(-nu*(y^tau)))^mu)^2)
                   d2ldvdt
                 }, 
                 
                 G.dev.incr = function(y, mu, sigma, nu, tau, ...) -2*dKW(y, mu, sigma, nu, tau, log = TRUE), 
                 rqres = expression(rqres(pfun = "pKW", type = "Continuous",  y = y, mu = mu, sigma = sigma, nu = nu, tau = tau)), 
                 
                 mu.initial = expression( mu <-  rep(0.5, length(y)) ), 
                 sigma.initial = expression( sigma <- rep(0.5, length(y)) ), 
                 nu.initial = expression( nu <- rep(0.5, length(y)) ), 
                 tau.initial = expression( tau <- rep(0.5, length(y)) ),
                 
                 mu.valid = function(mu) all(mu >  0), 
                 sigma.valid = function(sigma) all(sigma >  0), 
                 nu.valid = function(nu) all(nu > 0), 
                 tau.valid = function(tau) all(tau > 0),
                 
                 y.valid = function(y) all(y > 0)
  ), 
  class = c("gamlss.family", "family"))
}
#' @export
#' @rdname KW
dKW<-function(x,mu,sigma,nu,tau, log = FALSE){
  if (any(x<0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(mu<=0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma<=0)) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu<=0)) 
    stop(paste("nu must be positive", "\n", ""))
  if (any(tau<=0)) 
    stop(paste("tau must be positive", "\n", ""))
  
  term <- -mu*(x^sigma)
  loglik<- log(nu*tau*mu*sigma) + (sigma-1)*log(x) + term + 
    (nu-1)*log(1-exp(term)) + (tau-1)*log(1-(1-exp(term))^nu)
  
  if (log == FALSE) 
    density<- exp(loglik)
  else 
    density <- loglik
  return(density) 
}
#' @export
#' @rdname KW
pKW <- function(q,mu,sigma,nu,tau, lower.tail=TRUE, log.p = FALSE){
  if (any(q<0)) 
    stop(paste("q must be positive", "\n", ""))
  if (any(mu<=0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma<=0)) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu<=0)) 
    stop(paste("nu must be positive", "\n", ""))
  if (any(tau<=0)) 
    stop(paste("tau must be positive", "\n", ""))
  
  cdf <- 1-(1-(1-exp(-mu*(q^sigma)))^nu)^tau
  if (lower.tail == TRUE) 
    cdf <- cdf
  else cdf <- 1 - cdf
  if (log.p == FALSE) 
    cdf <- cdf
  else cdf <- log(cdf)
  cdf
  
}
#' @export
#' @rdname KW
qKW <- function(p,mu,sigma,nu,tau, lower.tail = TRUE, log.p = FALSE){
  
  if (any(mu<=0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma<=0)) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu<=0)) 
    stop(paste("nu must be positive", "\n", ""))
  if (any(tau<=0)) 
    stop(paste("tau must be positive", "\n", ""))
  
  if (log.p == TRUE) 
    p <- exp(p)
  else p <- p
  if (lower.tail == TRUE) 
    p <- p
  else p <- 1 - p
  if (any(p < 0) | any(p > 1)) 
    stop(paste("p must be between 0 and 1", "\n", ""))
  
  q <- ((-1/mu)*(log(1-(1-(1-p)^(1/tau))^(1/nu))))^(1/sigma)
  q
}
#' @export
#' @rdname KW
rKW <- function(n,mu,sigma,nu,tau){
  if(any(n<=0))
    stop(paste("n must be positive","\n",""))
  if (any(mu<=0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma<=0)) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu<=0)) 
    stop(paste("nu must be positive", "\n", ""))
  if (any(tau<=0)) 
    stop(paste("tau must be positive", "\n", ""))
  
  n <- ceiling(n)
  p <- runif(n)
  r <- qKW(p,mu,sigma,nu,tau)
  r
}
#' @export
#' @rdname KW
hKW<-function(x,mu,sigma,nu,tau){
  if (any(x<0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(mu<=0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma<=0)) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu<=0)) 
    stop(paste("nu must be positive", "\n", ""))
  if (any(tau<=0)) 
    stop(paste("tau must be positive", "\n", ""))
  
  h <- dKW(x,mu,sigma,nu,tau, log = FALSE)/pKW(q=x,mu,sigma,nu,tau, lower.tail=FALSE, log.p = FALSE)
  h
}




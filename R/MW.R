#' @name MW
#' 
#' @title 
#' The Modified Weibull Distribution
#' 
#' @description 
#' Density, distribution function, quantile function, 
#' random generation and hazard function for the modified weibull distribution with
#' parameters \code{mu}, \code{sigma} and \code{nu}.
#' 
#' @param x,q	vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. 
#' @param mu shape parameter one.    
#' @param sigma parameter two.
#' @param nu scale parameter three.        
#' @param log,log.p	logical; if TRUE, probabilities p are given as log(p).	
#' @param lower.tail logical; if TRUE (default), probabilities are 
#' P[X <= x], otherwise, P[X > x].
#' 
#' @details 
#' The modified weibull distribution with parameters \code{mu}, \code{sigma}
#' and \code{nu} has density given by
#' 
#' f(x) = mu*(sigma+nu*x)*x^(sigma-1)*exp(nu*x)*exp(-mu*x^(sigma)*exp(nu*x))
#' 
#' for x > 0.
#' 
#' @return 
#' \code{dMW} gives the density, \code{pMW} gives the distribution 
#' function, \code{qMW} gives the quantile function, \code{rMW}
#' generates random deviates and \code{hMW} gives the hazard function.
#' 
#' @export
#' @examples  
#' ## The probability density function 
#' curve(dMW(x, mu = 2, sigma = 1.5, nu = 0.2), from=0, to=2, ylim=c(0,1.5), col="red", las=1, ylab="The probability density function")
#' 
#' ## The cumulative distribution and the Reliability function
#' par(mfrow = c(1, 2))
#' curve(pMW(x, mu = 2, sigma = 1.5, nu = 0.2), from = 0, to = 2, ylim = c(0, 1), col = "red", las = 1, ylab = "The cumulative distribution function")
#' curve(pMW(x, mu = 2, sigma = 1.5, nu = 0.2, lower.tail = FALSE), from = 0, to = 2, ylim = c(0, 1), col = "red", las = 1, ylab = "The Reliability function")
#' 
#' ## The quantile function
#' p <- seq(from = 0, to = 0.998, length.out = 100)
#' plot(x = qMW(p = p, mu = 2, sigma = 1.5, nu = 0.2), y = p, xlab = "Quantile", las = 1, ylab = "Probability")
#' curve(pMW(x, mu = 2, sigma = 1.5, nu = 0.2), from = 0, add = TRUE, col = "red")
#' 
#' ## The random function
#' hist(rMW(n = 1000, mu = 2, sigma = 1.5, nu = 0.2), freq = FALSE, , ylim=c(0,1.5),xlab = "x", las = 1, main = "")
#' curve(dMW(x, mu = 2, sigma = 1.5, nu = 0.2), from = 0, , ylim=c(0,1.5), add = T, col = "red")
#' 
#' ## The Hazard function
#' curve(hMW(x, mu = 2, sigma = 1.5, nu = 0.2), from = 0, to = 1.5, ylim = c(0, 5), col = "red", las = 1, ylab = "The Hazard function")
MW <- function (mu.link = "log", sigma.link = "log", nu.link = "log") 
{
  mstats <- checklink("mu.link",    "Modified Weibull", substitute(mu.link),    c("log", "own"))
  dstats <- checklink("sigma.link", "Modified Weibull", substitute(sigma.link), c("log", "own"))
  vstats <- checklink("nu.link",    "Modified Weibull", substitute(nu.link),    c("log", "own"))
  
  structure(list(family = c("MW", "Modified Weibull"), 
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
                   dldm <- (1/mu)-y^sigma * exp(nu*y)  
                   dldm  
                 },
                 
                 d2ldm2 = function(y, mu, sigma, nu) {
                   d2ldm2<- -(1/(mu^2))
                   d2ldm2
                 }, 
                 
                 dldd = function(y, mu, sigma, nu) {
                   exp1<-mu*(y^sigma)*exp(nu*y)
                   dldd<-(1/(sigma+nu*y))+log(y)-log(y)*exp1
                   dldd          
                 },
                 
                 d2ldd2 = function(y, mu, sigma, nu) {
                   exp1<-mu*(y^sigma)*exp(nu*y)
                   dexp1dd<- y^sigma * log(y)*exp(nu*y)
                   d2ldd2<- -(1/(sigma+nu*y)^2) - log(y)*dexp1dd
                   d2ldd2
                 }, 
                 
                 dldv = function(y, mu, sigma, nu) {
                   exp1<-mu*(y^sigma)*exp(nu*y)
                   dldv<-y*((1/(sigma+nu*y))+1-exp1)
                   dldv
                 }, 
                 
                 d2ldv2 = function(y, mu, sigma, nu) {
                   exp1<-mu*(y^sigma)*exp(nu*y)
                   dexp1dv<- y*exp1
                   d2ldv2<- -y*((y/(sigma+nu*y)^2)+dexp1dv)
                   d2ldv2
                 }, 
                 
                 d2ldmdd = function(y, mu, sigma, nu) {
                   d2ldmdd<- -(y^sigma*log(y)*exp(nu*y))
                   d2ldmdd
                 }, 
                 
                 d2ldmdv = function(y, mu, sigma, nu) {
                   d2ldmdv<- -(y^(sigma+1)*exp(nu*y))
                   d2ldmdv
                 }, 
                 
                 d2ldddv = function(y, mu, sigma, nu) {
                   exp1<-mu*(y^sigma)*exp(nu*y)
                   dexp1dv<-  y*exp1
                   d2ldddv<- -(y/(sigma+nu*y)^2)-log(y)*dexp1dv
                   d2ldddv
                 }, 
                 
                 G.dev.incr = function(y, mu, sigma, nu,  ...) -2*dMW(y, mu, sigma, nu,  log = TRUE), 
                 rqres = expression(rqres(pfun = "pMW", type = "Continuous",  y = y, mu = mu, sigma = sigma, nu = nu)), 
                 
                 mu.initial = expression( mu <-  rep(1, length(y)) ), 
                 sigma.initial = expression( sigma <- rep(1, length(y)) ), 
                 nu.initial = expression( nu <- rep(1, length(y)) ), 
                 
                 mu.valid = function(mu) all(mu >  0), 
                 sigma.valid = function(sigma) all(sigma >=  0), 
                 nu.valid = function(nu) all(nu >= 0), 
                 
                 y.valid = function(y) all(y > 0)
  ), 
  class = c("gamlss.family", "family"))
}
#' @export
#' @rdname MW

dMW<-function(x,mu,sigma,nu, log = FALSE){
  if (any(x<0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(mu<=0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma<0)) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu<0)) 
    stop(paste("nu must be positive", "\n", ""))
  
  loglik<- log(mu) + log(sigma + nu*x) + (sigma-1)*log(x) +
    nu*x - mu*(x^sigma)*exp(nu*x)
  
  if (log == FALSE) 
    density<- exp(loglik)
  else 
    density <- loglik
  return(density)
}

#' @export
#' @rdname MW
pMW <- function(q,mu,sigma,nu, lower.tail=TRUE, log.p = FALSE){
  if (any(q<0)) 
    stop(paste("q must be positive", "\n", ""))
  if (any(mu<=0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma<0)) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu<0)) 
    stop(paste("nu must be positive", "\n", ""))
  
  cdf <- 1- exp(-mu*(q^sigma)*exp(nu*q))
  if (lower.tail == TRUE) 
    cdf <- cdf
  else cdf <- 1 - cdf
  if (log.p == FALSE) 
    cdf <- cdf
  else cdf <- log(cdf)
  cdf
}

#' @export
#' @rdname MW
qMW <- function(p, mu,sigma,nu,  lower.tail = TRUE, log.p = FALSE){
  if (any(mu<=0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma<0)) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu<0)) 
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
    1- exp(-mu*(x^sigma)*exp(nu*x))
  }
  
  fda1 <- function(x, mu,sigma,nu, p) {fda(x, mu,sigma,nu) - p}
  
  r_de_la_funcion <- function(mu,sigma,nu, p) {
    uniroot(fda1, interval=c(0,99999), mu,sigma,nu, p)$root
  }
  
  r_de_la_funcion <- Vectorize(r_de_la_funcion)
  q <- r_de_la_funcion(mu,sigma,nu, p)
  q
  
}

#' @export
#' @rdname MW
rMW <- function(n, mu,sigma,nu){
  if (any(mu<=0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma<0)) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu<0)) 
    stop(paste("nu must be positive", "\n", ""))
  
  n <- ceiling(n)
  p <- runif(n)
  r <- qMW(p,mu,sigma,nu) 
  r
}

#' @export
#' @rdname MW
hMW<-function(x,mu,sigma,nu){
  if (any(x<0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(mu<=0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma<0)) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu<0)) 
    stop(paste("nu must be positive", "\n", ""))
  
  h <- dMW(x,mu,sigma,nu, log = FALSE)/pMW(q=x,mu,sigma,nu, lower.tail=FALSE, log.p = FALSE)
  h  
}

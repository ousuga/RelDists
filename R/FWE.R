#' @name FWE
#' 
#' @title 
#' The Flexible Weibull Extension Distribution
#' 
#' @description 
#' Density, distribution function, quantile function, 
#' random generation  and hazard function for the flexible weibull extension distribution with
#' parameters \code{mu} and  \code{sigma}.
#' 
#' @param x,q	vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. 
#' @param mu parameter.    
#' @param sigma parameter.
#' @param log,log.p	logical; if TRUE, probabilities p are given as log(p).	
#' @param lower.tail logical; if TRUE (default), probabilities are 
#' P[X <= x], otherwise, P[X > x].
#' 
#' @details 
#' The flexible weibull extension with parameters \code{mu} and \code{sigma}
#' has density given by
#' 
#' f(x) = (mu + (sigma/x^2))*exp(mu*x - sigma/x)*exp(-exp(mu*x-sigma/x))
#' 
#' for x>0.
#' 
#' @return 
#' \code{dFWE} gives the density, \code{pFWE} gives the distribution 
#' function, \code{qFWE} gives the quantile function, \code{rFWE}
#' generates random deviates and \code{hFWE} gives the hazard function.
#' 
#' @export
#' @examples  
#' ## The probability density function
#' curve(dFWE(x, mu = 0.75, sigma = 0.5), from = 0, to = 3, ylim = c(0, 1.7), col = "red", las = 1, ylab = "The probability density function")
#' 
#' ## The cumulative distribution and the Reliability function
#' par(mfrow = c(1, 2))
#' curve(pFWE(x, mu = 0.75, sigma = 0.5), from = 0, to = 3, col = "red", las = 1, ylab = "The cumulative distribution function")
#' curve(pFWE(x, mu = 0.75, sigma = 0.5, lower.tail = FALSE), from = 0, to = 3, col = "red", las = 1, ylab = "The Reliability function")
#' 
#' ## The quantile function
#' p <- seq(from = 0, to = 0.99999, length.out = 100)
#' plot(x = qFWE(p,mu = 0.75, sigma = 0.5), y = p, xlab = "Quantile", las = 1, ylab = "Probability")
#' curve(pFWE(x, mu = 0.75, sigma = 0.5), from = 0, add = TRUE, col = "red")
#'  
#' ## The random function
#' hist(rFWE(n = 1000, mu = 2, sigma = 0.5), freq = FALSE, xlab = "x", ylim = c(0, 2), las = 1, main = "")
#' curve(dFWE(x, mu = 2, sigma = 0.5),  from = 0, to = 3, add = TRUE, , col = "red")
#' 
#' ## The Hazard function
#' curve(hFWE(x, mu = 0.75, sigma = 0.5), from = 0, to = 2, ylim = c(0, 2.5), col = "red", ylab = "The Hazard function", las = 1)
#' 
FWE <- function (mu.link="log", sigma.link="log") 
{
  mstats <- checklink("mu.link", "Flexible Weibull Extension", substitute(mu.link), c("log", "identity"))
  dstats <- checklink("sigma.link", "Flexible Weibull Extension", substitute(sigma.link), c("log", "identity"))
  
  structure(list(family = c("FEW", "Flexible Weibull Extension"),
                 parameters = list(mu=TRUE, sigma=TRUE), 
                 nopar = 2, 
                 type = "Continuous",
                 
                 mu.link = as.character(substitute(mu.link)), 
                 sigma.link = as.character(substitute(sigma.link)), 
                 
                 mu.linkfun = mstats$linkfun, 
                 sigma.linkfun = dstats$linkfun, 
                 
                 mu.linkinv = mstats$linkinv, 
                 sigma.linkinv = dstats$linkinv,
                 
                 mu.dr = mstats$mu.eta, 
                 sigma.dr = dstats$mu.eta,
                 
                 dldm = function(y,mu,sigma) (y-y*exp(mu*y-sigma/y)+y^2/(mu*y^2+sigma)),
                 
                 d2ldm2 = function(y,mu,sigma) {
                   dldm = function(y,mu,sigma) (y-y*exp(mu*y-sigma/y)+y^2/(mu*y^2+sigma))
                   ans <- dldm(y,mu,sigma)
                   ans <- -ans^2
                 },
                 
                 dldd = function(y,mu,sigma) (1/(mu*y^2+sigma)-1/y+exp(mu*y-sigma/y)/y),
                 
                 d2ldd2 = function(y,mu,sigma) {
                   dldd = function(y,mu,sigma) (1/(mu*y^2+sigma)-1/y+exp(mu*y-sigma/y)/y)
                   ans <- dldd(y,mu,sigma)
                   ans <- -ans^2
                 },
                 
                 d2ldmdd = function(y,mu,sigma) -(-y^2/(mu*y^2+sigma)^2+exp(mu*y-sigma/y))^2,
                 
                 G.dev.incr  = function(y,mu,sigma,...) -2*dFWE(y, mu, sigma, log=TRUE), 
                 rqres = expression(rqres(pfun="pFWE", type="Continuous", y=y, mu=mu, sigma=sigma)),
                 
                 mu.initial = expression( mu <-  rep(0.5, length(y)) ),     
                 sigma.initial = expression( sigma <- rep(0.5, length(y)) ), 
                 
                 mu.valid = function(mu) all(mu > 0) , 
                 sigma.valid = function(sigma)  all(sigma > 0), 
                 
                 y.valid = function(y)  all(y > 0)
  ),
  class = c("gamlss.family","family"))
}
#' @export
#' @rdname FWE
dFWE<-function(x,mu,sigma,log = FALSE){
  if (any(x<0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(mu<=0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma<=0)) 
    stop(paste("sigma must be positive", "\n", ""))
  
  loglik<- log(mu + (sigma/x^2)) + (mu*x) - (sigma/x) - 
    exp(mu*x - (sigma/x))
  
  if (log == FALSE) 
    density<- exp(loglik)
  else 
    density <- loglik
  return(density)
}
#' @export
#' @rdname FWE
pFWE <- function(q,mu,sigma, lower.tail=TRUE, log.p = FALSE){
  if (any(q<0)) 
    stop(paste("q must be positive", "\n", ""))
  if (any(mu<=0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma<=0)) 
    stop(paste("sigma must be positive", "\n", ""))
  
  cdf <- 1- exp(-exp(mu*q - sigma/q))
  
  if (lower.tail == TRUE) 
    cdf <- cdf
  else cdf <- 1 - cdf
  if (log.p == FALSE) 
    cdf <- cdf
  else cdf <- log(cdf)
  cdf
  
}
#' @export
#' @rdname FWE
qFWE <- function(p, mu, sigma, lower.tail = TRUE, log.p = FALSE) {
  if (any(mu<=0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma<=0)) 
    stop(paste("sigma must be positive", "\n", ""))
  
  if (log.p == TRUE) 
    p <- exp(p)
  else p <- p
  if (lower.tail == TRUE) 
    p <- p
  else p <- 1 - p
  if (any(p < 0) | any(p > 1)) 
    stop(paste("p must be between 0 and 1", "\n", ""))
  
  fda <- function(x,mu, sigma){
    1- exp(-exp(mu*x - sigma/x))
  }
  
  fda1 <- function(x, mu, sigma, p) {fda(x, mu, sigma) - p}
  
  r_de_la_funcion <- function(mu, sigma, p) {
    uniroot(fda1, interval=c(0,1e+06), mu, sigma, p)$root
  }
  
  r_de_la_funcion <- Vectorize(r_de_la_funcion)
  q <- r_de_la_funcion(mu, sigma, p)
  q
  
}

#' @export
#' @rdname FWE
rFWE <- function(n,mu,sigma){
  if (any(mu<=0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma<=0)) 
    stop(paste("sigma must be positive", "\n", ""))  
  
  n <- ceiling(n)
  p <- runif(n)
  r <- qFWE(p, mu,sigma)
  r
}
#' @export
#' @rdname FWE
hFWE<-function(x,mu,sigma){
  if (any(x<0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(mu<=0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma<=0)) 
    stop(paste("sigma must be positive", "\n", ""))
  
  h <- dFWE(x,mu,sigma, log = FALSE)/pFWE(q=x,mu,sigma, lower.tail=FALSE, log.p = FALSE)
  h
}



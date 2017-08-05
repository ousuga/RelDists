#' @name LW
#' 
#' @title 
#' The Log-Weibull Distribution
#' 
#' @description 
#' Density, distribution function, quantile function, 
#' random generation and hazard function for the log-weibull distribution with
#' parameters \code{mu} and \code{sigma}.
#' 
#' @param x,q	vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. 
#' @param mu parameter one.
#' @param sigma parameter two.
#' @param log,log.p	logical; if TRUE, probabilities p are given as log(p).	
#' @param lower.tail logical; if TRUE (default), probabilities are 
#' P[X <= x], otherwise, P[X > x].
#'  
#' @details 
#' The log-weibull distribution with parameters \code{mu} and 
#' \code{sigma} has density given by
#' 
#' f(x)= (1/sigma)*exp((x-mu)/sigma)*exp(-exp((x-mu)/sigma))
#' 
#' for -Inf < x < Inf.
#' 
#' @return 
#' \code{dLW} gives the density, \code{pLW} gives the distribution 
#' function, \code{qLW} gives the quantile function, \code{rLW}
#' generates random deviates and \code{hLW} gives the hazard function.
#' 
#' @export
#' @examples  
#' ## The probability density function 
#' curve(dLW(x, mu = 0, sigma = 1), from = -20, to = 10, ylim = c(0, 0.4), col = "red", las = 1, ylab = "The probability density function")
#' 
#' ## The cumulative distribution and the Reliability function
#' par(mfrow = c(1, 2))
#' curve(pLW(x, mu = 0, sigma = 1), from = 0, to = 10, ylim = c(0, 1), col = "red", las = 1, ylab = "The cumulative distribution function")
#' curve(pLW(x, mu = 0, sigma = 1, lower.tail = FALSE), from = 0, to = 10, ylim = c(0, 1), col = "red", las = 1, ylab = "The Reliability function")
#' 
#' ## The quantile function
#' p <- seq(from = 0, to = 0.998, length.out = 100)
#' plot(x=qLW(p, mu = 0, sigma = 1), y = p, xlab = "Quantile", las = 1, ylab = "Probability")
#' curve(pLW(x, mu = 0, sigma = 1), from = -6, add = TRUE, col = "red")
#' 
#' ## The random function
#' hist(rLW(10000, mu = 0, sigma = 1), freq = FALSE,  ylim = c(0, 0.4), xlab = "x", las = 1, main = "")
#' curve(dLW(x, mu = 0, sigma = 1),  from = -20, to = 10,  ylim = c(0, 0.4), add = TRUE, col = "red") 
#' 
#' ## The Hazard function
#' curve(hLW(x, mu = 0, sigma = 1), from = -20, to = 0, ylim = c(0, 0.3), col = "red", ylab = "The hazard function", las = 1)
#'
LW <- function (mu.link="identity", sigma.link="log") 
{
  mstats <- checklink("mu.link", "Log-Weibull", substitute(mu.link), c("identity", "own"))
  dstats <- checklink("sigma.link", "Log-Weibull", substitute(sigma.link), c("log", "own"))
  
  structure(list(family = c("LW", "Log-Weibull"),
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
                 
                 dldm = function(y,mu,sigma) (-(1/sigma) + (1/sigma)*exp((y-mu)/sigma)),
                 d2ldm2 = function(y,mu,sigma) -(-(1/sigma^2)*exp((y-mu)/sigma))^2,
                 dldd = function(y,mu,sigma) (-(1/sigma) -(y/sigma^2) + (y/sigma^2)*exp((y-mu)/sigma)),
                 d2ldd2 = function(y,mu,sigma) -((1/sigma^2) + (2*y/sigma^3) -(2*y/sigma^3)*exp((y-mu)/sigma) - ((y/sigma^2)^2)*exp((y-mu)/sigma) )^2,
                 d2ldmdd = function(y,mu,sigma) -((1/sigma^2) - (1/sigma^2)*exp((y-mu)/sigma) -(y/sigma^3)*exp((y-mu)/sigma))^2,
                 
                 G.dev.incr  = function(y,mu,sigma,...) -2*dLW(y, mu, sigma, log=TRUE), 
                 rqres = expression(rqres(pfun="pLW", type="Continuous", y=y, mu=mu, sigma=sigma)),
                 
                 mu.initial = expression( mu <-  rep(0.5, length(y)) ),     
                 sigma.initial = expression( sigma <- rep(0.5, length(y)) ), 
                 
                 mu.valid = function(mu) TRUE , 
                 sigma.valid = function(sigma)  all(sigma > 0), 
                 
                 y.valid = function(y)  TRUE),
            
            class = c("gamlss.family","family"))
}
#' @export
#' @rdname LW
dLW<-function(x,mu,sigma, log = FALSE){
  if (any(sigma<=0)) 
    stop(paste("sigma must be positive", "\n", ""))
  
  loglik<- -log(sigma) + (x-mu)/sigma - exp((x-mu)/sigma)
  
  if (log == FALSE) 
    density<- exp(loglik)
  else 
    density <- loglik
  return(density)
}

#' @export
#' @rdname LW
pLW <- function(q,mu,sigma, lower.tail=TRUE, log.p = FALSE){
  if (any(sigma<=0)) 
    stop(paste("sigma must be positive", "\n", ""))
  
  cdf <- 1-exp(-exp((q-mu)/sigma))
  
  if (lower.tail == TRUE) 
    cdf <- cdf
  else cdf <- 1 - cdf
  if (log.p == FALSE) 
    cdf <- cdf
  else cdf <- log(cdf)
  cdf
}
#' @export
#' @rdname LW
qLW <- function(p,mu,sigma, lower.tail = TRUE, log.p = FALSE){
  
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
  
  q <- sigma*(log(-log(1-p)))+mu  
  q
}
#' @export
#' @rdname LW
rLW <- function(n,mu,sigma){
  if(any(n<=0))
    stop(paste("n must be positive","\n",""))
  if (any(sigma<=0)) 
    stop(paste("sigma must be positive", "\n", ""))
  
  n <- ceiling(n)
  p <- runif(n)
  r <- qLW(p,mu,sigma)
  r
}
#' @export
#' @rdname LW
hLW<-function(x,mu,sigma){
  if (any(sigma<=0)) 
    stop(paste("sigma must be positive", "\n", ""))
  
  h <- dLW(x,mu,sigma, log = FALSE)/pLW(q=x,mu,sigma, lower.tail=FALSE, log.p = FALSE)
  h
}

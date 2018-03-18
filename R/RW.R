#' @name RW
#' 
#' @title 
#' The Reflected Weibull Distribution
#' 
#' @description 
#' Density, distribution function, quantile function, 
#' random generation and hazard function for the reflected weibull distribution with
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
#' The reflected weibull distribution with parameters \code{mu} and
#' \code{sigma} has density given by
#' 
#' f(x) = mu*sigma*(-x)^(sigma-1)*exp(-mu*(-x)^sigma)
#' 
#' for - inf < x < 0.  
#' 
#' @return 
#' \code{dRW} gives the density, \code{pRW} gives the distribution 
#' function, \code{qRW} gives the quantile function, \code{rRW}
#' generates random deviates and \code{hRW} gives the hazard function.
#' 
#' @export
#' @examples  
#' ## The probability density function
#' curve(dRW(x, mu = 1, sigma = 1), from = -5, to = 0, ylim = c(0, 1), col = "red", las = 1, ylab = "The probability density function")
#' 
#' ## The cumulative distribution and the Reliability function
#' par(mfrow = c(1, 2))
#' curve(pRW(x, mu = 1, sigma = 1), from = -5, to = 0, ylim = c(0, 1), col = "red", las = 1, ylab ="The cumulative distribution function")
#' curve(pRW(x, mu = 1, sigma = 1, lower.tail = FALSE), from = -5, to = 0, ylim = c(0, 1), col = "red", las = 1, ylab = "The Reliability function")
#' 
#' ## The quantile function
#' p <- seq(from = 0, to = 0.99999, length.out = 100)
#' plot(x=qRW(p=p,mu = 1, sigma = 1), y=p, xlab="Quantile", las=1, ylab="Probability")
#' curve(pRW(x, mu = 1, sigma = 1), from = -5, add = TRUE, col = "red")
#' 
#' ## The random function
#' hist(rRW(n = 10000, mu = 1, sigma = 1), freq = FALSE,xlab = "x", las = 1, main = "")
#' curve(dRW(x, mu = 1, sigma = 1),  from = -5, to = 0, add = TRUE, col = "red")
#' 
#' ## The Hazard function
#' curve(hRW(x, mu = 1, sigma = 1), from = -5, to = 0, ylim = c(0, 1), col = "red", ylab = "The hazard function", las = 1)
RW <- function (mu.link="log", sigma.link="log") 
{
  mstats <- checklink("mu.link"   , "Reflected Weibull", substitute(mu.link), c("log", "own"))
  dstats <- checklink("sigma.link", "Reflected Weibull", substitute(sigma.link), c("log", "own"))
  
  structure(list(family = c("RW", "Reflected Weibull"),
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
                 
                 dldm = function(y,mu,sigma) ((1/mu) -((-y)^sigma) ),
                 d2ldm2 = function(mu) (-1/mu^2),
                 dldd = function(y,mu,sigma) ((1/sigma)+log(-y)-mu*(-((-y)^sigma)*log(-y))),
                 d2ldd2 = function(y,mu,sigma) {
                   dldd = function(y,mu,sigma) ((-1/(sigma^2))-mu*((-y)^sigma)*(log(-y))^2)
                   ans <- dldd(y,mu,sigma)
                   ans <- -ans^2
                 },
                 d2ldmdd = function(y,mu,sigma) -(((-y)^sigma)*log(-y))^2,
                 
                 G.dev.incr  = function(y,mu,sigma,...) -2*dRW(y, mu, sigma, log=TRUE), 
                 rqres = expression(rqres(pfun="pRW", type="Continuous", y=y, mu=mu, sigma=sigma)),
                 
                 mu.initial = expression( mu <-  rep(0.5, length(y)) ),     
                 sigma.initial = expression( sigma <- rep(0.5, length(y)) ), 
                 
                 mu.valid = function(mu) all(mu > 0) , 
                 sigma.valid = function(sigma)  all(sigma > 0), 
                 
                 y.valid = function(y)  all(y < 0)
  ),
  class = c("gamlss.family","family"))
}
#' @export
#' @rdname RW

dRW<-function(x,mu,sigma, log=FALSE){
  if (any(x>0)) 
    stop(paste("x must be negative", "\n", ""))
  if (any(mu <= 0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma<=0)) 
    stop(paste("sigma must be positive", "\n", ""))
  
  loglik<- log(mu) + log(sigma) + (sigma-1)*log(-x) -
    mu*((-x)^sigma)
  
  if (log == FALSE) 
    density<- exp(loglik)
  else 
    density <- loglik
  return(density)
}

#' @export
#' @rdname RW
pRW <- function(q,mu,sigma, lower.tail=TRUE, log.p = FALSE){
  # if (any(q<0)) 
  #  stop(paste("q must be positive", "\n", ""))
  if (any(mu <= 0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma<=0)) 
    stop(paste("sigma must be positive", "\n", ""))
  
  cdf <- exp(-mu*(-q)^sigma)
  
  if (lower.tail == TRUE) 
    cdf <- cdf
  else cdf <- 1 - cdf
  if (log.p == FALSE) 
    cdf <- cdf
  else cdf <- log(cdf)
  cdf
}

#' @export
#' @rdname RW
qRW <- function(p,mu,sigma, lower.tail = TRUE, log.p = FALSE){
  if (any(mu <= 0 )) 
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
  
  q <- -{((-1/mu)*log(p))^(1/sigma)}
  q
}

#' @export
#' @rdname RW
rRW <- function(n,mu,sigma){
  if(any(n<=0))
    stop(paste("n must be positive","\n",""))
  if (any(mu<=0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma<=0)) 
    stop(paste("sigma must be positive", "\n", ""))
  
  n <- ceiling(n)
  p <- runif(n)
  r <- qRW(p,mu,sigma)
  r
}
#' @export
#' @rdname RW
hRW<-function(x,mu,sigma){
  if (any(x>0)) 
    stop(paste("x must be negative", "\n", ""))
  if (any(mu <= 0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma<=0)) 
    stop(paste("sigma must be positive", "\n", ""))
  
  h <- dRW(x,mu,sigma, log = FALSE)/pRW(q=x,mu,sigma, lower.tail=FALSE, log.p = FALSE)
  h
}

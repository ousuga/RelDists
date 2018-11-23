#' @name IW
#' 
#' @title 
#' The Inverse Weibull Distribution
#' 
#' @description 
#' Density, distribution function, quantile function, 
#' random generation and hazard function for the inverse weibull distribution with
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
#' The inverse weibull distribution with parameters \code{mu} and
#' \code{sigma} has density given by
#' 
#' f(x) = mu*sigma*x^(-sigma-1)*exp(-mu*(x^-sigma)),
#' 
#' for x > 0.
#' 
#' @return 
#' \code{dIW} gives the density, \code{pIW} gives the distribution 
#' function, \code{qIW} gives the quantile function, \code{rIW}
#' generates random deviates and \code{hIW} gives the hazard function.
#' 
#' @export
#' @examples  
#' ## The probability density function 
#' curve(dIW(x, mu = 5, sigma = 2.5), from = 0, to = 10, ylim = c(0, 0.55), col = "red", las = 1, ylab = "The probability density function")
#' 
#' ## The cumulative distribution and the Reliability function
#' par(mfrow = c(1, 2))
#' curve(pIW(x, mu = 5, sigma = 2.5), from = 0, to = 10, ylim = c(0, 1), col = "red", las = 1, ylab = "The cumulative distribution function")
#' curve(pIW(x, mu = 5, sigma = 2.5, lower.tail = FALSE), from = 0, to = 10, ylim = c(0, 1), col = "red", las = 1, ylab = "The Reliability function")
#' 
#' ## The quantile function
#' p <- seq(from = 0, to = 0.998, length.out = 100)
#' plot(x = qIW(p, mu = 5, sigma = 2.5), y = p, xlab = "Quantile", las = 1, ylab = "Probability")
#' curve(pIW(x, mu = 5, sigma = 2.5),  from = 0, add = TRUE, col = "red")
#' 
#' ## The random function
#' hist(rIW(n = 1000, mu = 5, sigma = 2.5), freq = FALSE, xlab = "x", ylim = c(0, 0.55), las = 1, main = "")
#' curve(dIW(x, mu = 5, sigma = 2.5),  from = 0, to = 10, add = TRUE, ylim = c(0, 0.55), col = "red")
#' 
#' ## The Hazard function
#' curve(hIW(x, mu = 5, sigma = 2.5), from = 0, to = 15, ylim = c(0, 1), col = "red", las = 1, ylab = "The Hazard function")
IW <- function (mu.link="log", sigma.link="log") 
{
  mstats <- checklink("mu.link", "Inverse Weibull", substitute(mu.link), c("log", "own"))
  dstats <- checklink("sigma.link", "Inverse Weibull", substitute(sigma.link), c("log", "own"))
  
  structure(list(family = c("IW", "Inverse Weibull"),
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
                 
                 dldm = function(y,mu,sigma) (1/mu)-y^(-sigma),
                 d2ldm2 = function(mu) -1/(mu^2),
                 dldd = function(y,mu,sigma) (1/sigma)-log(y) +mu*(y^(-sigma))*log(y),
                 d2ldd2 = function(y,mu,sigma) {
                   dldd = function(y,mu,sigma) -(1/(sigma^2))-mu*((log(y))^2)*(y^(-sigma))
                   ans <- dldd(y,mu,sigma)
                   ans <- -ans^2
                 },
                 
                 d2ldmdd = function(y,sigma) -(log(y)*(y^(-sigma)))^2,
                 
                 G.dev.incr  = function(y,mu,sigma,...) -2*dIW(y, mu, sigma, log=TRUE), 
                 rqres = expression(rqres(pfun="pIW", type="Continuous", y=y, mu=mu, sigma=sigma)),
                 
                 mu.initial = expression( mu <-  rep(0.5, length(y)) ),     
                 sigma.initial = expression( sigma <- rep(0.5, length(y)) ),
                 
                 mu.valid = function(mu) all(mu > 0) , 
                 sigma.valid = function(sigma)  all(sigma > 0), 
                 
                 y.valid = function(y)  all(y > 0)
  ),
  class = c("gamlss.family","family"))
}
#' @export
#' @rdname IW
dIW<-function(x,mu,sigma, log = FALSE){
  if (any(x<0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(mu <= 0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma<=0)) 
    stop(paste("sigma must be positive", "\n", ""))
  
  loglik<- log(mu) + log(sigma) - (sigma+1)*log(x) - 
    mu*(x^-sigma)
  
  if (log == FALSE) 
    density<- exp(loglik)
  else 
    density <- loglik
  return(density)  
}

#' @export
#' @rdname IW
pIW <- function(q,mu,sigma, lower.tail=TRUE, log.p = FALSE){
  if (any(q<0)) 
    stop(paste("q must be positive", "\n", ""))
  if (any(mu<=0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma<=0)) 
    stop(paste("sigma must be positive", "\n", ""))
  
  cdf <- exp((-mu)*(q^(-sigma)))
  
  if (lower.tail == TRUE) 
    cdf <- cdf
  else cdf <- 1 - cdf
  if (log.p == FALSE) 
    cdf <- cdf
  else cdf <- log(cdf)
  cdf
}

#' @export
#' @rdname IW
qIW <- function(p,mu,sigma, lower.tail = TRUE, log.p = FALSE){
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
  
  q <- ((-1/mu)*log(p))^(-1/sigma)
  q
}

#' @export
#' @rdname IW
rIW <- function(n,mu,sigma){
  if(any(n<=0))
    stop(paste("n must be positive","\n",""))
  if (any(mu<=0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma<=0)) 
    stop(paste("sigma must be positive", "\n", ""))
  
  n <- ceiling(n)
  p <- runif(n)
  r <- qIW(p, mu,sigma)
  r
}

#' @export
#' @rdname IW
hIW<-function(x,mu,sigma){
  if (any(x<0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(mu <= 0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma<=0)) 
    stop(paste("sigma must be positive", "\n", ""))
  
  h <- dIW(x,mu,sigma, log = FALSE)/pIW(q=x,mu,sigma, lower.tail=FALSE, log.p = FALSE)
  h
}



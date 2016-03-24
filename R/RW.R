#' @name RW
#' 
#' @title 
#' The Reflected Weibull Distribution
#' 
#' @description 
#' Density, distribution function, quantile function, 
#' random generation and hazard function for the reflected weibull distribution with
#' parameters \code{alpha} and \code{theta}.
#' 
#' @param x,q	vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. 
#' @param alpha parameter one.
#' @param theta parameter two.   
#' @param log,log.p	logical; if TRUE, probabilities p are given as log(p).	
#' @param lower.tail logical; if TRUE (default), probabilities are 
#' P[X <= x], otherwise, P[X > x].
#'
#' @details 
#' The reflected weibull distribution with parameters \code{alpha} and
#' \code{theta} has density given by
#' 
#' f(x) = alpha*theta*(-x)^(theta-1)*exp(-alpha*(-x)^alpha)
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
#' curve(dRW(x, alpha = 1, theta = 1), from = -5, to = 0, ylim = c(0, 1), col = "red", las = 1, ylab = "The probability density function")
#' 
#' ## The cumulative distribution and the Reliability function
#' par(mfrow = c(1, 2))
#' curve(pRW(x, alpha = 1, theta = 1), from = -5, to = 0, ylim = c(0, 1), col = "red", las = 1, ylab ="The cumulative distribution function")
#' curve(pRW(x, alpha = 1, theta = 1, lower.tail = FALSE), from = -5, to = 0, ylim = c(0, 1), col = "red", las = 1, ylab = "The Reliability function")
#' 
#' ## The quantile function
#' p <- seq(from = 0, to = 0.99999, length.out = 100)
#' plot(x=qRW(p=p,alpha = 1, theta = 1), y=p, xlab="Quantile", las=1, ylab="Probability")
#' curve(pRW(x, alpha = 1, theta = 1), from = -5, add = TRUE, col = "red")
#' 
#' ## The random function
#' hist(rRW(n = 10000, alpha = 1, theta = 1), freq = FALSE,xlab = "x", las = 1, main = "")
#' curve(dRW(x, alpha = 1, theta = 1),  from = -5, to = 0, add = TRUE, col = "red")
#' 
#' ## The Hazard function
#' curve(hRW(x, alpha = 1, theta = 1), from = -5, to = 0, ylim = c(0, 1), col = "red", ylab = "The hazard function", las = 1)


dRW<-function(x,alpha,theta, log=FALSE){
  if (any(x>0)) 
    stop(paste("x must be negative", "\n", ""))
  if (any(alpha <= 0 )) 
    stop(paste("alpha must be positive", "\n", ""))
  if (any(theta<=0)) 
    stop(paste("theta must be positive", "\n", ""))
  
  loglik<- log(alpha) + log(theta) + (theta-1)*log(-x) -
    alpha*((-x)^theta)
  
  if (log == FALSE) 
    density<- exp(loglik)
  else 
    density <- loglik
  return(density)
}

#' @export
#' @rdname RW
pRW <- function(q,alpha,theta, lower.tail=TRUE, log.p = FALSE){
  # if (any(q<0)) 
  #  stop(paste("q must be positive", "\n", ""))
  if (any(alpha <= 0 )) 
    stop(paste("alpha must be positive", "\n", ""))
  if (any(theta<=0)) 
    stop(paste("theta must be positive", "\n", ""))
  
  cdf <- exp(-alpha*(-q)^theta)
  
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
qRW <- function(p,alpha,theta, lower.tail = TRUE, log.p = FALSE){
  if (any(alpha <= 0 )) 
    stop(paste("alpha must be positive", "\n", ""))
  if (any(theta<=0)) 
    stop(paste("theta must be positive", "\n", ""))
  
  if (log.p == TRUE) 
    p <- exp(p)
  else p <- p
  if (lower.tail == TRUE) 
    p <- p
  else p <- 1 - p
  if (any(p < 0) | any(p > 1)) 
    stop(paste("p must be between 0 and 1", "\n", ""))
  
  q <- -{((-1/alpha)*log(p))^(1/theta)}
  q
}

#' @export
#' @rdname RW
rRW <- function(n,alpha,theta){
  if(any(n<=0))
    stop(paste("n must be positive","\n",""))
  if (any(alpha<=0 )) 
    stop(paste("alpha must be positive", "\n", ""))
  if (any(theta<=0)) 
    stop(paste("theta must be positive", "\n", ""))
  
  n <- ceiling(n)
  p <- runif(n)
  r <- qRW(p,alpha,theta)
  r
}
#' @export
#' @rdname RW
hRW<-function(x,alpha,theta){
  if (any(x>0)) 
    stop(paste("x must be negative", "\n", ""))
  if (any(alpha <= 0 )) 
    stop(paste("alpha must be positive", "\n", ""))
  if (any(theta<=0)) 
    stop(paste("theta must be positive", "\n", ""))
  
  h <- dRW(x,alpha,theta, log = FALSE)/pRW(q=x,alpha,theta, lower.tail=FALSE, log.p = FALSE)
  h
}

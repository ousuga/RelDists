#' @name MWEx
#' 
#' @title 
#' The Modified Weibull Extension Distribution 
#' 
#' @description 
#' Density, distribution function, quantile function, 
#' random generation and hazard function for the modified weibull extension distribution with
#' parameters \code{alpha}, \code{beta} and \code{lambda}.
#' 
#' @param x,q	vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. 
#' @param alpha parameter one.    
#' @param beta parameter two.
#' @param lambda parameter three.        
#' @param log,log.p	logical; if TRUE, probabilities p are given as log(p).	
#' @param lower.tail logical; if TRUE (default), probabilities are 
#' P[X <= x], otherwise, P[X > x].
#' 
#' @details 
#' The modified weibull extension distribution with parameters \code{alpha}, \code{beta}
#' and \code{lambda} has density given by
#' 
#' f(x) = lambda*beta*((x/alpha)^(beta-1))*exp((x/alpha)^beta + lambda*alpha*(1-exp((x/alpha)^beta)))
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
#' curve(dMWEx(x, alpha = 1/0.5, beta = 3, lambda = 2), from = 0, to = 2.5, ylim = c(0, 1.5), col = "red", las = 1, ylab = "The probability density function")
#' 
#' ## The cumulative distribution and the Reliability function
#' par(mfrow = c(1, 2))
#' curve(pMWEx(x, alpha = 1/0.5, beta = 3, lambda = 2), from = 0, to = 2.5, ylim = c(0, 1), col = "red", las = 1, ylab = "The cumulative distribution function")
#' curve(pMWEx(x, alpha = 1/0.5, beta = 3, lambda = 2,  lower.tail = FALSE), from = 0, to = 2.5, ylim = c(0, 1), col = "red", las = 1, ylab = "The Reliability function")
#' 
#' ## The quantile function
#' p <- seq(from = 0, to = 0.998, length.out = 100)
#' plot(x = qMWEx(p, alpha = 1/0.5, beta = 3, lambda = 2), y = p, xlab = "Quantile", las = 1, ylab = "Probability")
#' curve(pMWEx(x, alpha = 1/0.5, beta = 3, lambda = 2), from = 0, add = TRUE, col = "red")
#' 
#' ## The random function
#' hist(rMWEx(n = 10000, alpha = 1/0.5, beta = 3, lambda = 2), freq = FALSE, ylim = c(0, 1.5), xlab = "x", las = 1, main = "")
#' curve(dMWEx(x, alpha = 1/0.5, beta = 3, lambda = 2),  from = 0, ylim = c(0, 2.5), add = T, col = "red")
#' 
#' ## The Hazard function
#' curve(hMWEx(x, alpha = 1/0.5, beta = 3, lambda = 2), from = 0, to = 1.7, ylim = c(0, 12), col = "red", ylab = "The hazard function", las = 1)


dMWEx<-function(x,alpha,beta,lambda, log = FALSE){
  if (any(x<0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(alpha<=0 )) 
    stop(paste("alpha must be positive", "\n", ""))
  if (any(beta<=0)) 
    stop(paste("beta must be positive", "\n", ""))
  if (any(lambda<=0)) 
    stop(paste("lambda must be positive", "\n", ""))
  
  
  loglik<-log(lambda*beta) + (beta-1)*log(x/alpha) +
    (x/alpha)^beta + lambda*alpha*(1-exp((x/alpha)^beta))
  
  if (log == FALSE) 
    density<- exp(loglik)
  else 
    density <- loglik
  return(density)
}

#' @export
#' @rdname MWEx

pMWEx<- function(q,alpha,beta,lambda, lower.tail=TRUE, log.p = FALSE){
  if (any(q<0)) 
    stop(paste("q must be positive", "\n", ""))
  if (any(alpha<=0 )) 
    stop(paste("alpha must be positive", "\n", ""))
  if (any(beta<=0)) 
    stop(paste("beta must be positive", "\n", ""))
  if (any(lambda<=0)) 
    stop(paste("lambda must be positive", "\n", ""))
  
  cdf <- 1- exp(lambda*alpha*(1-exp((q/alpha)^beta)))
  
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
qMWEx <- function(p,alpha,beta,lambda, lower.tail = TRUE, log.p = FALSE){
  
  if (any(alpha<=0 )) 
    stop(paste("alpha must be positive", "\n", ""))
  if (any(beta<=0)) 
    stop(paste("beta must be positive", "\n", ""))
  if (any(lambda<=0)) 
    stop(paste("lambda must be positive", "\n", ""))
  
  if (log.p == TRUE) 
    p <- exp(p)
  else p <- p
  if (lower.tail == TRUE) 
    p <- p
  else p <- 1 - p
  if (any(p < 0) | any(p > 1)) 
    stop(paste("p must be between 0 and 1", "\n", ""))
  
  q <- alpha*(log(1-((1/(lambda*alpha))*log(1-p))))^(1/beta)
  q
}


#' @export
#' @rdname MWEx
rMWEx <- function(n,alpha,beta,lambda){
  if(any(n<=0))
    stop(paste("n must be positive","\n",""))
  if (any(alpha<=0 )) 
    stop(paste("alpha must be positive", "\n", ""))
  if (any(beta<=0)) 
    stop(paste("beta must be positive", "\n", ""))
  if (any(lambda<=0)) 
    stop(paste("lambda must be positive", "\n", ""))
  n <- ceiling(n)
  p <- runif(n)
  r <- qMWEx(p, alpha,beta,lambda)
  r
}

#' @export
#' @rdname MWEx
hMWEx<-function(x,alpha,beta,lambda){
  if (any(x<0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(alpha<=0 )) 
    stop(paste("alpha must be positive", "\n", ""))
  if (any(beta<=0)) 
    stop(paste("beta must be positive", "\n", ""))
  if (any(lambda<=0)) 
    stop(paste("lambda must be positive", "\n", ""))
  
  h <- dMWEx(x,alpha,beta,lambda, log = FALSE)/pMWEx(q=x,alpha,beta,lambda, lower.tail=FALSE, log.p = FALSE)
  h
}

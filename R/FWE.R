#' @name FWE
#' 
#' @title 
#' The Flexible Weibull Extension Distribution
#' 
#' @description 
#' Density, distribution function, quantile function, 
#' random generation  and hazard function for the flexible weibull extension distribution with
#' parameters \code{alpha} and  \code{beta}.
#' 
#' @param x,q	vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. 
#' @param alpha parameter.    
#' @param beta parameter.
#' @param log,log.p	logical; if TRUE, probabilities p are given as log(p).	
#' @param lower.tail logical; if TRUE (default), probabilities are 
#' P[X <= x], otherwise, P[X > x].
#' 
#' @details 
#' The flexible weibull extension with parameters \code{alpha} and \code{beta}
#' has density given by
#' 
#' f(x) = (alpha + (beta/x^2))*exp(alpha*x - beta/x)*exp(-exp(alpha*x-beta/x))
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
#' curve(dFWE(x, alpha = 0.75, beta = 0.5), from = 0, to = 3, ylim = c(0, 1.7), col = "red", las = 1, ylab = "The probability density function")
#' 
#' ## The cumulative distribution and the Reliability function
#' par(mfrow = c(1, 2))
#' curve(pFWE(x, alpha = 0.75, beta = 0.5), from = 0, to = 3, col = "red", las = 1, ylab = "The cumulative distribution function")
#' curve(pFWE(x, alpha = 0.75, beta = 0.5, lower.tail = FALSE), from = 0, to = 3, col = "red", las = 1, ylab = "The Reliability function")
#' 
#' ## The quantile function
#' p <- seq(from = 0, to = 0.99999, length.out = 100)
#' plot(x = qFWE(p,alpha = 0.75, beta = 0.5), y = p, xlab = "Quantile", las = 1, ylab = "Probability")
#' curve(pFWE(x, alpha = 0.75, beta = 0.5), from = 0, add = TRUE, col = "red")
#'  
#' ## The random function
#' hist(rFWE(n = 1000, alpha = 2, beta = 0.5), freq = FALSE, xlab = "x", ylim = c(0, 2), las = 1, main = "")
#' curve(dFWE(x, alpha = 2, beta = 0.5),  from = 0, to = 3, add = TRUE, , col = "red")
#' 
#' ## The Hazard function
#' curve(hFWE(x, alpha = 0.75, beta = 0.5), from = 0, to = 2, ylim = c(0, 2.5), col = "red", ylab = "The Hazard function", las = 1)
#' 
dFWE<-function(x,alpha,beta,log = FALSE){
  if (any(x<0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(alpha<=0 )) 
    stop(paste("alpha must be positive", "\n", ""))
  if (any(beta<=0)) 
    stop(paste("beta must be positive", "\n", ""))
  
  loglik<- log(alpha + (beta/x^2)) + (alpha*x) - (beta/x) - 
    exp(alpha*x - (beta/x))
  
  if (log == FALSE) 
    density<- exp(loglik)
  else 
    density <- loglik
  return(density)
}
#' @export
#' @rdname FWE
pFWE <- function(q,alpha,beta, lower.tail=TRUE, log.p = FALSE){
  if (any(q<0)) 
    stop(paste("q must be positive", "\n", ""))
  if (any(alpha<=0 )) 
    stop(paste("alpha must be positive", "\n", ""))
  if (any(beta<=0)) 
    stop(paste("beta must be positive", "\n", ""))
  
  cdf <- 1- exp(-exp(alpha*q - beta/q))
  
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
qFWE <- function(p, alpha, beta, lower.tail = TRUE, log.p = FALSE) {
  if (any(alpha<=0 )) 
    stop(paste("alpha must be positive", "\n", ""))
  if (any(beta<=0)) 
    stop(paste("beta must be positive", "\n", ""))
  
  if (log.p == TRUE) 
    p <- exp(p)
  else p <- p
  if (lower.tail == TRUE) 
    p <- p
  else p <- 1 - p
  if (any(p < 0) | any(p > 1)) 
    stop(paste("p must be between 0 and 1", "\n", ""))
  
  fda <- function(x,alpha, beta){
    1- exp(-exp(alpha*x - beta/x))
  }
  
  fda1 <- function(x, alpha, beta, p) {fda(x, alpha, beta) - p}
  
  r_de_la_funcion <- function(alpha, beta, p) {
    uniroot(fda1, interval=c(0,1e+06), alpha, beta, p)$root
  }
  
  r_de_la_funcion <- Vectorize(r_de_la_funcion)
  q <- r_de_la_funcion(alpha, beta, p)
  q
  
}

#' @export
#' @rdname FWE
rFWE <- function(n,alpha,beta){
  if (any(alpha<=0 )) 
    stop(paste("alpha must be positive", "\n", ""))
  if (any(beta<=0)) 
    stop(paste("beta must be positive", "\n", ""))  
  
  n <- ceiling(n)
  p <- runif(n)
  r <- qFWE(p, alpha,beta)
  r
}
#' @export
#' @rdname FWE
hFWE<-function(x,alpha,beta){
  if (any(x<0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(alpha<=0 )) 
    stop(paste("alpha must be positive", "\n", ""))
  if (any(beta<=0)) 
    stop(paste("beta must be positive", "\n", ""))
  
  h <- dFWE(x,alpha,beta, log = FALSE)/pFWE(q=x,alpha,beta, lower.tail=FALSE, log.p = FALSE)
  h
}



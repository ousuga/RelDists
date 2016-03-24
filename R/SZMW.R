#' @name SZMW
#' 
#' @title 
#' The Sarhan and Zaindins Modified Weibull Distribution
#' 
#' @description 
#' Density, distribution function, quantile function, 
#' random generation and hazard function for Sarhan and Zaindins modified weibull distribution with
#' parameters \code{alpha}, \code{beta} and \code{gamma}.
#' 
#' @param x,q	vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. 
#' @param alpha parameter one.    
#' @param beta parameter two.
#' @param gamma parameter three.
#' @param log,log.p	logical; if TRUE, probabilities p are given as log(p).	
#' @param lower.tail logical; if TRUE (default), probabilities are 
#' P[X <= x], otherwise, P[X > x].
#' 
#' @details 
#' The Sarhan and Zaindins modified weibull with parameters \code{alpha}, 
#' \code{beta} and \code{gamma} has density given by
#' 
#' f(x)=(alpha+beta*gamma*x^(gamma-1))*exp(-alpha*x-beta*x^gamma)
#' 
#' for x>0.
#' 
#' @return 
#' \code{dSZMW} gives the density, \code{pSZMW} gives the distribution 
#' function, \code{qSZMW} gives the quantile function, \code{rSZMW}
#' generates random deviates and \code{hSZMW} gives the hazard function.
#'  
#' @export
#' @examples 
#' 
#' ## The probability density function
#' curve(dSZMW(x, alpha = 2, beta = 1.5, gamma = 0.2), from = 0, to = 2, ylim = c(0, 1.7), col = "red", las = 1, ylab = "The probability density function")
#' 
#' ## The cumulative distribution and the Reliability function
#' par(mfrow = c(1, 2))
#' curve(pSZMW(x, alpha = 2, beta = 1.5, gamma = 0.2), from = 0, to = 2, ylim = c(0, 1), col = "red", las = 1, ylab = "The cumulative distribution function")
#' curve(pSZMW(x, alpha = 2, beta = 1.5, gamma = 0.2, lower.tail = FALSE), from = 0, to = 2, ylim = c(0, 1), col = "red", las = 1, ylab = "The Reliability function")
#' 
#' ## The quantile function
#' p <- seq(from = 0, to = 0.99999, length.out = 100)
#' plot(x = qSZMW(p = p, alpha = 2, beta = 1.5, gamma = 0.2), y = p, xlab = "Quantile", las = 1, ylab = "Probability")
#' curve(pSZMW(x, alpha = 2, beta = 1.5, gamma = 0.2), from = 0, add = TRUE, col = "red")
#' 
#' ## The random function
#' hist(rSZMW(n = 1000, alpha = 2, beta = 1.5, gamma = 0.2), freq = FALSE, xlab = "x", las = 1, main = "")
#' curve(dSZMW(x, alpha = 2, beta = 1.5, gamma = 0.2),  from = 0, add = TRUE, col = "red")
#' 
#' ## The Hazard function
#' curve(hSZMW(x, alpha = 2, beta = 1.5, gamma = 0.2), from = 0, to = 3, ylim = c(0, 8), col = "red", ylab = "The hazard function", las = 1)
#' 
dSZMW<-function(x,alpha,beta,gamma, log = FALSE){
  if (any(x<0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(alpha<=0 )) 
    stop(paste("alpha must be positive", "\n", ""))
  if (any(beta<=0)) 
    stop(paste("beta must be positive", "\n", ""))
  if (any(gamma<=0)) 
    stop(paste("gamma must be positive", "\n", ""))
  
  loglik<- log(alpha + beta*gamma*x^(gamma-1)) - alpha*x - beta*x^gamma
  
  if (log == FALSE) 
    density<- exp(loglik)
  else 
    density <- loglik
  return(density)
}
#' @export
#' @rdname SZMW
pSZMW <- function(q,alpha,beta,gamma, lower.tail=TRUE, log.p = FALSE){
  if (any(q<0)) 
    stop(paste("q must be positive", "\n", ""))
  if (any(alpha<=0 )) 
    stop(paste("alpha must be positive", "\n", ""))
  if (any(beta<=0)) 
    stop(paste("beta must be positive", "\n", ""))
  if (any(gamma<=0)) 
    stop(paste("gamma must be positive", "\n", ""))
  
  cdf <- 1- exp(-alpha*q -beta*(q^gamma))
  if (lower.tail == TRUE) 
    cdf <- cdf
  else cdf <- 1 - cdf
  if (log.p == FALSE) 
    cdf <- cdf
  else cdf <- log(cdf)
  cdf
  
}
#' @export
#' @rdname SZMW
qSZMW <- function(p, alpha,beta,gamma, lower.tail = TRUE, log.p = FALSE) {
  if (any(alpha<=0 )) 
    stop(paste("alpha must be positive", "\n", ""))
  if (any(beta<=0)) 
    stop(paste("beta must be positive", "\n", ""))
  if (any(gamma<=0)) 
    stop(paste("gamma must be positive", "\n", ""))
  
  if (log.p == TRUE) 
    p <- exp(p)
  else p <- p
  if (lower.tail == TRUE) 
    p <- p
  else p <- 1 - p
  if (any(p < 0) | any(p > 1)) 
    stop(paste("p must be between 0 and 1", "\n", ""))
  
  fda <- function(x,alpha,beta,gamma){
    1- exp(-alpha*x - beta*(x^gamma))
  }
  
  fda1 <- function(x, alpha,beta,gamma, p) {fda(x, alpha,beta,gamma) - p}
  
  r_de_la_funcion <- function(alpha,beta,gamma, p) {
    uniroot(fda1, interval=c(0,1e+06), alpha,beta,gamma, p)$root
  }
  
  r_de_la_funcion <- Vectorize(r_de_la_funcion)
  q <- r_de_la_funcion(alpha,beta,gamma, p)
  q
  
}
#' @export
#' @rdname SZMW
rSZMW<- function(n,alpha,beta,gamma){
  if (any(alpha<=0 )) 
    stop(paste("alpha must be positive", "\n", ""))
  if (any(beta<=0)) 
    stop(paste("beta must be positive", "\n", ""))
  if (any(gamma<=0)) 
    stop(paste("gamma must be positive", "\n", ""))
  
  n <- ceiling(n)
  p <- runif(n)
  r <- qSZMW(p, alpha,beta,gamma)
  r
}
#' @export
#' @rdname SZMW
hSZMW<-function(x,alpha,beta,gamma){
  if (any(x<0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(alpha<=0 )) 
    stop(paste("alpha must be positive", "\n", ""))
  if (any(beta<=0)) 
    stop(paste("beta must be positive", "\n", ""))
  if (any(gamma<=0)) 
    stop(paste("gamma must be positive", "\n", ""))
  
  h <- dSZMW(x,alpha,beta,gamma, log = FALSE)/pSZMW(q=x,alpha,beta,gamma, lower.tail=FALSE, log.p = FALSE)
  h  
}



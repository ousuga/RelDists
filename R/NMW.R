#' @name NMW
#' 
#' @title 
#' The Almaki and Yuan's modified Weibull distribution
#' 
#' @description 
#' Density, distribution function, quantile function, 
#' random generation and hazard function for the Almaki and Yuan's modified weibull distribution with
#' parameters \code{alpha}, \code{beta}, \code{theta}, \code{gamma} and \code{lambda}.
#' 
#' @param x,q	vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. 
#' @param alpha parameter one.    
#' @param beta parameter two.
#' @param theta parameter three.
#' @param gamma parameter four.
#' @param lambda parameter five.
#' @param log,log.p	logical; if TRUE, probabilities p are given as log(p).	
#' @param lower.tail logical; if TRUE (default), probabilities are 
#' P[X <= x], otherwise, P[X > x].
#' 
#' @details 
#' The Almaki and Yuans modified weibull with parameters \code{alpha}, 
#' \code{beta}, \code{theta}, \code{gamma} and \code{lambda} has density given by
#' 
#' f(x)=(alpha*theta*x^(theta-1)+beta*(gamma+ lambda*x)*(x^(gamma-1)*exp(lambda*x)))*(exp((-alpha*x^theta)-(beta*x^gamma*exp(lambda*x)))
#' 
#' for x>0.
#' 
#' @return 
#' \code{dNMW} gives the density, \code{pNMW} gives the distribution 
#' function, \code{qNMW} gives the quantile function, \code{rNMW}
#' generates random deviates and \code{hNMW} gives the hazard function.
#' 
#' @export
#' @examples 
#' ## The probability density function 
#' curve(dNMW(x, alpha = 1.15, beta = 0.15, theta = 0.75, gamma = 5, lambda = 2), from = 0, to = 1.4, ylim = c(0, 3), col = "red", las = 1, ylab = "The probability density function")
#' 
#' ## The cumulative distribution and the Reliability function
#' par(mfrow = c(1, 2))
#' curve(pNMW(x, alpha = 1.15, beta = 0.15, theta = 0.75, gamma = 5, lambda = 2), from = 0, to = 1.4, col = "red", las = 1, ylab = "The cumulative distribution function")
#' curve(pNMW(x, alpha = 1.15, beta = 0.15, theta = 0.75, gamma = 5, lambda = 2, lower.tail = FALSE), from = 0, to = 1.4, col = "red", las = 1, ylab = "The Reliability function")
#' 
#' ## The quantile function
#' p <- seq(from = 0, to = 0.998, length.out = 100)
#' plot(x=qNMW(p, alpha = 1.15, beta = 0.15, theta = 0.75, gamma = 5, lambda = 2), y = p, xlab = "Quantile", las = 1, ylab = "Probability")
#' curve(pNMW(x, alpha = 1.15, beta = 0.15, theta = 0.75, gamma = 5, lambda = 2), from = 0, add = TRUE, col = "red")
#' 
#' ## The random function
#' hist(rNMW(n = 1000, alpha = 1.15, beta = 0.15, theta = 0.75, gamma = 5, lambda = 2), freq = FALSE, ylim = c(0, 3), xlab = "x", las = 1, main = "")
#' curve(dNMW(x, alpha = 1.15, beta = 0.15, theta = 0.75, gamma = 5, lambda = 2),  from = 0, ylim = c(0, 3), add = T, col = "red")
#' 
#' ## The Hazard function
#' curve(hNMW(x, alpha = 1.15, beta = 0.15, theta = 0.75, gamma = 5, lambda = 2), from = 0, to = 1.5, ylim = c(0, 8), col = "red", ylab = "The hazard function", las = 1)
#' 
dNMW<-function(x,alpha,beta,theta,gamma,lambda, log = FALSE){
  if (any(x<0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(alpha<=0 )) 
    stop(paste("alpha must be positive", "\n", ""))
  if (any(beta<=0)) 
    stop(paste("beta must be positive", "\n", ""))
  if (any(theta<=0)) 
    stop(paste("theta must be positive", "\n", ""))
  if (any(gamma<=0)) 
    stop(paste("gamma must be positive", "\n", ""))
  if (any(lambda<0)) 
    stop(paste("lambda must be positive", "\n", ""))
  
  loglik<- log(alpha*theta*(x^(theta-1))  + exp(lambda*x)*beta*(gamma+lambda*x)*x^(gamma-1)) -
    alpha*(x^theta) - beta*(x^gamma)*exp(lambda*x) 
  
  if (log == FALSE) 
    density<- exp(loglik)
  else 
    density <- loglik
  return(density)
}
#' @export
#' @rdname NMW
pNMW <- function(q,alpha,beta,theta,gamma,lambda, lower.tail=TRUE, log.p = FALSE){
  if (any(q<0)) 
    stop(paste("q must be positive", "\n", ""))
  if (any(alpha<=0 )) 
    stop(paste("alpha must be positive", "\n", ""))
  if (any(beta<=0)) 
    stop(paste("beta must be positive", "\n", ""))
  if (any(theta<=0)) 
    stop(paste("theta must be positive", "\n", ""))
  if (any(gamma<=0)) 
    stop(paste("gamma must be positive", "\n", ""))
  if (any(lambda<0)) 
    stop(paste("lambda must be positive", "\n", ""))
  
  cdf <- 1-exp(-alpha*(q^theta) -beta*(q^gamma)*exp(lambda*q))
  
  if (lower.tail == TRUE) 
    cdf <- cdf
  else cdf <- 1 - cdf
  if (log.p == FALSE) 
    cdf <- cdf
  else cdf <- log(cdf)
  cdf
}
#' @export
#' @rdname NMW
qNMW <- function(p,alpha,beta,theta,gamma,lambda, lower.tail = TRUE, log.p = FALSE){
  if (any(alpha<=0 )) 
    stop(paste("alpha must be positive", "\n", ""))
  if (any(beta<=0)) 
    stop(paste("beta must be positive", "\n", ""))
  if (any(theta<=0)) 
    stop(paste("theta must be positive", "\n", ""))
  if (any(gamma<=0)) 
    stop(paste("gamma must be positive", "\n", ""))
  if (any(lambda<0)) 
    stop(paste("lambda must be positive", "\n", ""))
  
  if (log.p == TRUE) 
    p <- exp(p)
  else p <- p
  if (lower.tail == TRUE) 
    p <- p
  else p <- 1 - p
  if (any(p < 0) | any(p > 1)) 
    stop(paste("p must be between 0 and 1", "\n", ""))
  
  
  fda <- function(x,alpha,beta,theta,gamma,lambda){
    1-exp(-alpha*(x^theta) -beta*(x^gamma)*exp(lambda*x))
  }
  
  fda1 <- function(x, alpha,beta,theta,gamma,lambda, p) {fda(x, alpha,beta,theta,gamma,lambda) - p}
  
  r_de_la_funcion <- function(alpha,beta,theta,gamma,lambda, p) {
    uniroot(fda1, interval=c(0,1e+06), alpha,beta,theta,gamma,lambda, p)$root
  }
  
  r_de_la_funcion <- Vectorize(r_de_la_funcion)
  q <- r_de_la_funcion(alpha,beta,theta,gamma,lambda, p)
  q
}
#' @export
#' @rdname NMW
rNMW <- function(n,alpha,beta,theta,gamma,lambda){
  if(any(n<=0))
    stop(paste("n must be positive","\n",""))
  if (any(alpha<=0 )) 
    stop(paste("alpha must be positive", "\n", ""))
  if (any(beta<=0)) 
    stop(paste("beta must be positive", "\n", ""))
  if (any(theta<=0)) 
    stop(paste("theta must be positive", "\n", ""))
  if (any(gamma<=0)) 
    stop(paste("gamma must be positive", "\n", ""))
  if (any(lambda<0)) 
    stop(paste("lambda must be positive", "\n", ""))
  
  n <- ceiling(n)
  p <- runif(n)
  r <- qNMW(p,alpha,beta,theta,gamma,lambda)
  r
}
#' @export
#' @rdname NMW
hNMW<-function(x,alpha,beta,theta,gamma,lambda, log = FALSE){
  if (any(x<0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(alpha<=0 )) 
    stop(paste("alpha must be positive", "\n", ""))
  if (any(beta<=0)) 
    stop(paste("beta must be positive", "\n", ""))
  if (any(theta<=0)) 
    stop(paste("theta must be positive", "\n", ""))
  if (any(gamma<=0)) 
    stop(paste("gamma must be positive", "\n", ""))
  if (any(lambda<0)) 
    stop(paste("lambda must be positive", "\n", ""))
  
  h <- dNMW(x,alpha,beta,theta,gamma,lambda, log = FALSE)/pNMW(q=x,alpha,beta,theta,gamma,lambda, lower.tail=FALSE, log.p = FALSE)
  h 
}


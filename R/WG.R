#' @name WG
#' 
#' 
#' @title 
#' The Weibull Geometric Distribution
#' 
#' @description 
#' Density, distribution function, quantile function, 
#' random generation and hazard function for the weibull geometric distribution with
#' parameters \code{alpha}, \code{beta} and \code{pi}.
#' 
#' @param x,q	vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. 
#' @param alpha shape parameter.
#' @param beta scale parameter.
#' @param pi parameter of geometric random variable.             
#' @param log,log.p	logical; if TRUE, probabilities p are given as log(p).	
#' @param lower.tail logical; if TRUE (default), probabilities are 
#' P[X <= x], otherwise, P[X > x].
#'  
#' @details 
#' The weibull geometric distribution with parameters \code{alpha},
#' \code{beta} and \code{pi} has density given by
#' 
#' f(x) = (alpha*(beta)^alpha*(1-pi)*x^(alpha-1)*exp(-(beta*x)^alpha))/(1- pi*exp(-(beta*x)^alpha))^2
#' 
#' for x > 0.
#' 
#' @return 
#' \code{dWG} gives the density, \code{pWG} gives the distribution 
#' function, \code{qWG} gives the quantile function, \code{rWG}
#' generates random deviates and \code{hWG} gives the hazard function.
#' 
#' @export
#' @examples  
#' ## The probability density function 
#' curve(dWG(x, alpha = 0.5, beta = 0.2, pi = 0.95), from = 0, to = 5, ylim = c(0, 1.5), col = "red", las = 1, ylab = "The probability density function")
#' 
#' ## The cumulative distribution and the Reliability function
#' par(mfrow = c(1, 2))
#' curve(pWG(x, alpha = 0.5, beta = 0.2, pi = 0.95), from = 0, to = 5, ylim = c(0, 1), col = "red", las = 1, ylab = "The cumulative distribution function")
#' curve(pWG(x, alpha = 0.5, beta = 0.2, pi = 0.95, lower.tail = FALSE), from = 0, to = 5, ylim = c(0, 1), col = "red", las = 1, ylab = "The Reliability function")
#' 
#' ## The quantile function
#' p <- seq(from = 0, to = 0.99999, length.out = 100)
#' plot(x = qWG(p = p, alpha = 0.5, beta = 0.2, pi = 0.95), y = p, xlab = "Quantile", las = 1, ylab = "Probability")
#' curve(pWG(x, alpha = 0.5, beta = 0.2, pi = 0.95), from = 0, add = TRUE, col = "red")
#' 
#' ## The random function
#' hist(rWG(1000, alpha = 0.5, beta = 0.2, pi = 0.95), freq = FALSE, xlab = "x", ylim = c(0, 0.2), las = 1, main = "")
#' curve(dWG(x, alpha = 0.5, beta = 0.2, pi = 0.95),  from = 0, add = TRUE, col = "red", ylim = c(0, 0.2))
#' 
#' ## The Hazard function
#' curve(hWG(x, alpha = 0.5, beta = 0.2, pi = 0.95), from = 0, to = 15, ylim = c(0, 2.2), col = "red", ylab = "The hazard function", las = 1)


dWG<-function(x,alpha,beta,pi, log = FALSE){
  if (any(x<0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(alpha<=0 )) 
    stop(paste("alpha must be positive", "\n", ""))
  if (any(beta<=0)) 
    stop(paste("beta must be positive", "\n", ""))
  if (any(pi <=0  | pi>=1  )) 
    stop(paste("pi must be between zero and one", "\n", ""))
  
  loglik<- log(alpha) + alpha*log(beta) + log(1-pi) + (alpha-1)*log(x) - 
    (beta*x)^alpha - 2*log(1-pi*exp(-(beta*x)^alpha))
  if (log == FALSE) 
    density<- exp(loglik)
  else 
    density <- loglik
  return(density)
}

#' @export
#' @rdname WG
pWG <- function(q,alpha,beta,pi, lower.tail=TRUE, log.p = FALSE){
  if (any(q<0)) 
    stop(paste("q must be positive", "\n", ""))
  if (any(alpha<=0 )) 
    stop(paste("alpha must be positive", "\n", ""))
  if (any(beta<=0)) 
    stop(paste("beta must be positive", "\n", ""))
  if (any(pi <=0  | pi>=1  )) 
    stop(paste("p must be between zero and one", "\n", ""))
  
  cdf <- (1-exp(-(beta*q)^alpha))/(1-pi*exp(-(beta*q)^alpha))
  
  if (lower.tail == TRUE) 
    cdf <- cdf
  else cdf <- 1 - cdf
  if (log.p == FALSE) 
    cdf <- cdf
  else cdf <- log(cdf)
  cdf
}

#' @export
#' @rdname WG
qWG <- function(p, alpha, beta, pi, lower.tail = TRUE, log.p = FALSE) {
  if (any(alpha<=0 )) 
    stop(paste("alpha must be positive", "\n", ""))
  if (any(beta<=0)) 
    stop(paste("beta must be positive", "\n", ""))
  if (any(pi <=0  | pi>=1  )) 
    stop(paste("pi must be between zero and one", "\n", ""))
  
  if (log.p == TRUE) 
    p <- exp(p)
  else p <- p
  if (lower.tail == TRUE) 
    p <- p
  else p <- 1 - p
  if (any(p < 0) | any(p > 1)) 
    stop(paste("p must be between 0 and 1", "\n", ""))
  
  
  fda <- function(x,alpha, beta,pi){
    (1- exp(-(beta*x)^alpha))/(1-(pi*exp(-(beta*x)^alpha)))
  }
  
  fda1 <- function(x, alpha, beta, pi, p) {fda(x, alpha, beta,pi) - p}
  
  r_de_la_funcion <- function(alpha, beta, pi,p) {
    uniroot(fda1, interval=c(0,1e+06), alpha, beta, pi,p)$root
  }
  
  r_de_la_funcion <- Vectorize(r_de_la_funcion)
  q <- r_de_la_funcion(alpha, beta, pi,p)
  q
  
}

#' @export
#' @rdname WG
rWG <- function(n,alpha,beta,pi){
  if (any(alpha<=0 )) 
    stop(paste("alpha must be positive", "\n", ""))
  if (any(beta<=0)) 
    stop(paste("beta must be positive", "\n", ""))  
  if (any(pi <=0  | pi>=1  )) 
    stop(paste("pi must be between zero and one", "\n", ""))
  
  n <- ceiling(n)
  p <- runif(n)
  r <- qWG(p, alpha,beta,pi)
  r
}
#' @export
#' @rdname WG
hWG<-function(x,alpha,beta,pi){
  if (any(x<0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(alpha<=0 )) 
    stop(paste("alpha must be positive", "\n", ""))
  if (any(beta<=0)) 
    stop(paste("beta must be positive", "\n", ""))  
  if (any(pi <=0  | pi>=1  )) 
    stop(paste("pi must be between zero and one", "\n", ""))
  
  h <- dWG(x,alpha,beta,pi, log = FALSE)/pWG(q=x,alpha,beta,pi, lower.tail=FALSE, log.p = FALSE)
  h  
}

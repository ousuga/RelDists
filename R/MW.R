#' @name MW
#' 
#' @title 
#' The Modified Weibull Distribution
#' 
#' @description 
#' Density, distribution function, quantile function, 
#' random generation and hazard function for the modified weibull distribution with
#' parameters \code{beta}, \code{gamma} and \code{lambda}.
#' 
#' @param x,q	vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. 
#' @param beta shape parameter one.    
#' @param gamma parameter two.
#' @param lambda scale parameter three.        
#' @param log,log.p	logical; if TRUE, probabilities p are given as log(p).	
#' @param lower.tail logical; if TRUE (default), probabilities are 
#' P[X <= x], otherwise, P[X > x].
#' 
#' @details 
#' The modified weibull distribution with parameters \code{beta}, \code{gamma}
#' and \code{lambda} has density given by
#' 
#' f(x) = beta*(gamma+lambda*x)*x^(gamma-1)*exp(lambda*x)*exp(-beta*x^(gamma)*exp(lambda*x))
#' 
#' for x > 0.
#' 
#' @return 
#' \code{dMW} gives the density, \code{pMW} gives the distribution 
#' function, \code{qMW} gives the quantile function, \code{rMW}
#' generates random deviates and \code{hMW} gives the hazard function.
#' 
#' @export
#' @examples  
#' ## The probability density function 
#' curve(dMW(x, beta = 2, gamma = 1.5, lambda = 0.2), from=0, to=2, ylim=c(0,1.5), col="red", las=1, ylab="The probability density function")
#' 
#' ## The cumulative distribution and the Reliability function
#' par(mfrow = c(1, 2))
#' curve(pMW(x, beta = 2, gamma = 1.5, lambda = 0.2), from = 0, to = 2, ylim = c(0, 1), col = "red", las = 1, ylab = "The cumulative distribution function")
#' curve(pMW(x, beta = 2, gamma = 1.5, lambda = 0.2, lower.tail = FALSE), from = 0, to = 2, ylim = c(0, 1), col = "red", las = 1, ylab = "The Reliability function")
#' 
#' ## The quantile function
#' p <- seq(from = 0, to = 0.998, length.out = 100)
#' plot(x = qMW(p = p, beta = 2, gamma = 1.5, lambda = 0.2), y = p, xlab = "Quantile", las = 1, ylab = "Probability")
#' curve(pMW(x, beta = 2, gamma = 1.5, lambda = 0.2), from = 0, add = TRUE, col = "red")
#' 
#' ## The random function
#' hist(rMW(n = 1000, beta = 2, gamma = 1.5, lambda = 0.2), freq = FALSE, , ylim=c(0,1.5),xlab = "x", las = 1, main = "")
#' curve(dMW(x, beta = 2, gamma = 1.5, lambda = 0.2), from = 0, , ylim=c(0,1.5), add = T, col = "red")
#' 
#' ## The Hazard function
#' curve(hMW(x, beta = 2, gamma = 1.5, lambda = 0.2), from = 0, to = 1.5, ylim = c(0, 5), col = "red", las = 1, ylab = "The Hazard function")


dMW<-function(x,beta,gamma,lambda, log = FALSE){
  if (any(x<0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(beta<=0 )) 
    stop(paste("beta must be positive", "\n", ""))
  if (any(gamma<0)) 
    stop(paste("gamma must be positive", "\n", ""))
  if (any(lambda<0)) 
    stop(paste("lambda must be positive", "\n", ""))
  
  loglik<- log(beta) + log(gamma + lambda*x) + (gamma-1)*log(x) +
    lambda*x - beta*(x^gamma)*exp(lambda*x)
  
  if (log == FALSE) 
    density<- exp(loglik)
  else 
    density <- loglik
  return(density)
}

#' @export
#' @rdname MW
pMW <- function(q,beta,gamma,lambda, lower.tail=TRUE, log.p = FALSE){
  if (any(q<0)) 
    stop(paste("q must be positive", "\n", ""))
  if (any(beta<=0 )) 
    stop(paste("beta must be positive", "\n", ""))
  if (any(gamma<0)) 
    stop(paste("gamma must be positive", "\n", ""))
  if (any(lambda<0)) 
    stop(paste("lambda must be positive", "\n", ""))
  
  cdf <- 1- exp(-beta*(q^gamma)*exp(lambda*q))
  if (lower.tail == TRUE) 
    cdf <- cdf
  else cdf <- 1 - cdf
  if (log.p == FALSE) 
    cdf <- cdf
  else cdf <- log(cdf)
  cdf
}

#' @export
#' @rdname MW
qMW <- function(p, beta,gamma,lambda,  lower.tail = TRUE, log.p = FALSE){
  if (any(beta<=0 )) 
    stop(paste("beta must be positive", "\n", ""))
  if (any(gamma<0)) 
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
  
  fda <- function(x,beta,gamma,lambda){
    1- exp(-beta*(x^gamma)*exp(lambda*x))
  }
  
  fda1 <- function(x, beta,gamma,lambda, p) {fda(x, beta,gamma,lambda) - p}
  
  r_de_la_funcion <- function(beta,gamma,lambda, p) {
    uniroot(fda1, interval=c(0,99999), beta,gamma,lambda, p)$root
  }
  
  r_de_la_funcion <- Vectorize(r_de_la_funcion)
  q <- r_de_la_funcion(beta,gamma,lambda, p)
  q
  
}

#' @export
#' @rdname MW
rMW <- function(n, beta,gamma,lambda){
  if (any(beta<=0 )) 
    stop(paste("beta must be positive", "\n", ""))
  if (any(gamma<0)) 
    stop(paste("gamma must be positive", "\n", ""))
  if (any(lambda<0)) 
    stop(paste("lambda must be positive", "\n", ""))
  
  n <- ceiling(n)
  p <- runif(n)
  r <- qMW(p,beta,gamma,lambda) 
  r
}

#' @export
#' @rdname MW
hMW<-function(x,beta,gamma,lambda){
  if (any(x<0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(beta<=0 )) 
    stop(paste("beta must be positive", "\n", ""))
  if (any(gamma<0)) 
    stop(paste("gamma must be positive", "\n", ""))
  if (any(lambda<0)) 
    stop(paste("lambda must be positive", "\n", ""))
  
  h <- dMW(x,beta,gamma,lambda, log = FALSE)/pMW(q=x,beta,gamma,lambda, lower.tail=FALSE, log.p = FALSE)
  h  
}

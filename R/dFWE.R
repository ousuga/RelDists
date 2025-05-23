#' The Flexible Weibull Extension distribution
#' 
#' @description
#' Density, distribution function, quantile function, 
#' random generation and hazard function for the Flexible Weibull Extension distribution with
#' parameters \code{mu} and \code{sigma}.
#' 
#' @param x,q	vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. 
#' @param mu parameter.    
#' @param sigma parameter.
#' @param log,log.p	logical; if TRUE, probabilities p are given as log(p).	
#' @param lower.tail logical; if TRUE (default), probabilities are 
#' P[X <= x], otherwise, P[X > x].
#' 
#' @seealso \link{FWE}
#' 
#' @details 
#' The Flexible Weibull extension with parameters \code{mu} and \code{sigma}
#' has density given by
#' 
#' \eqn{f(x) = (\mu + \sigma/x^2) \exp(\mu x - \sigma/x) \exp(-\exp(\mu x-\sigma/x))}
#' 
#' for x>0.
#' 
#' @return 
#' \code{dFWE} gives the density, \code{pFWE} gives the distribution 
#' function, \code{qFWE} gives the quantile function, \code{rFWE}
#' generates random deviates and \code{hFWE} gives the hazard function.
#' 
#' @example examples/examples_dFWE.R
#' 
#' @export
dFWE <- function(x, mu, sigma, log=FALSE){
  if (any(x < 0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(mu <= 0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  
  loglik <- log(mu + (sigma/x^2)) + (mu*x) - (sigma/x) - 
    exp(mu*x - (sigma/x))
  
  if (log == FALSE) 
    density <- exp(loglik)
  else density <- loglik
  return(density)
}
#' @export
#' @rdname dFWE
pFWE <- function(q, mu, sigma, lower.tail=TRUE, log.p=FALSE){
  if (any(q < 0)) 
    stop(paste("q must be positive", "\n", ""))
  if (any(mu <= 0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  
  cdf <- 1- exp(-exp(mu*q - sigma/q))
  
  if (lower.tail == TRUE) 
    cdf <- cdf
  else cdf <- 1 - cdf
  if (log.p == FALSE) 
    cdf <- cdf
  else cdf <- log(cdf)
  cdf
}
#' @importFrom stats uniroot
#' @export
#' @rdname dFWE
qFWE <- function(p, mu, sigma, lower.tail=TRUE, log.p=FALSE) {
  if (any(mu <= 0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  
  if (log.p == TRUE) 
    p <- exp(p)
  else p <- p
  if (lower.tail == TRUE) 
    p <- p
  else p <- 1 - p
  if (any(p < 0) | any(p > 1)) 
    stop(paste("p must be between 0 and 1", "\n", ""))
  q <- (log(-log(1-p)) + sqrt((log(-log(1-p)))^2 + 4*mu*sigma)) / (2*mu)
  q
}
#' @importFrom stats runif
#' @export
#' @rdname dFWE
rFWE <- function(n, mu, sigma){
  if (any(mu <= 0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  
  n <- ceiling(n)
  p <- runif(n)
  r <- qFWE(p, mu, sigma)
  r
}
#' @export
#' @rdname dFWE
hFWE <- function(x, mu, sigma){
  if (any(x < 0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(mu <= 0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  
  h <- (mu+sigma/x^2) * exp(mu*x-sigma/x)
  h
}
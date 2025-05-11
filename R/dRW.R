#' The Reflected Weibull distribution
#' 
#' @author Amylkar Urrea Montoya, \email{amylkar.urrea@@udea.edu.co}
#' 
#' @description 
#' Density, distribution function, quantile function, 
#' random generation and hazard function for the Reflected Weibull Distribution 
#' with parameters \code{mu} and \code{sigma}.
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
#' @seealso \link{RW}
#' 
#' @details 
#' The Reflected Weibull Distribution with parameters \code{mu} 
#' and \code{sigma} has density given by
#' 
#' \eqn{f(x) = \mu \sigma (-x) ^{\sigma - 1} e^{-\mu(-x)^\sigma},}
#' 
#' for \eqn{x < 0}.
#' 
#' @return 
#' \code{dRW} gives the density, \code{pRW} gives the distribution 
#' function, \code{qRW} gives the quantile function, \code{rRW}
#' generates random deviates and \code{hRW} gives the hazard function.
#'
#' @example examples/examples_dRW.R
#'
#' @references
#' Almalki, S. J., & Nadarajah, S. (2014). Modifications of the 
#' Weibull distribution: A review. Reliability Engineering & 
#' System Safety, 124, 32-55.
#' 
#' Cohen, A. C. (1973). The reflected Weibull distribution. 
#' Technometrics, 15(4), 867-873.
#'
#' @export
dRW <- function(x, mu, sigma, log=FALSE){
  if (any(x >= 0)) 
    stop(paste("x must be negative", "\n", ""))
  if (any(mu <= 0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  
  loglik<- log(mu) + log(sigma) + (sigma-1)*log(-x) -
    mu*((-x)^sigma)
  
  if (log == FALSE) 
    density <- exp(loglik)
  else 
    density <- loglik
  return(density)
}
#' @export
#' @rdname dRW
pRW <- function(q, mu, sigma,
                lower.tail=TRUE, log.p=FALSE){
  if (any(q >= 0)) 
    stop(paste("q must be negative", "\n", ""))
  if (any(mu <= 0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
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
#' @rdname dRW
qRW <- function(p, mu, sigma,
                lower.tail=TRUE, log.p=FALSE){
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
  
  q <- -{((-1/mu)*log(p))^(1/sigma)}
  q
}
#' @importFrom stats runif
#' @export
#' @rdname dRW
rRW <- function(n, mu, sigma){
  if(any(n <= 0))
    stop(paste("n must be positive","\n",""))
  if (any(mu <= 0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  
  n <- ceiling(n)
  p <- runif(n)
  r <- qRW(p, mu, sigma)
  r
}
#' @export
#' @rdname dRW
hRW <- function(x, mu, sigma){
  if (any(x >= 0)) 
    stop(paste("x must be negative", "\n", ""))
  if (any(mu <= 0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  
  h <- dRW(x, mu, sigma, log=FALSE) /
    pRW(q=x, mu, sigma, lower.tail=FALSE, log.p=FALSE)
  h
}


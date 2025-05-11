#' The Inverse Weibull distribution
#'
#' @author Johan David Marin Benjumea, \email{johand.marin@udea.edu.co}
#' 
#' @description 
#' Density, distribution function, quantile function, 
#' random generation and hazard function for the inverse weibull distribution with
#' parameters \code{mu} and \code{sigma}.
#' 
#' @param x,q	vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. 
#' @param mu scale parameter.
#' @param sigma shape parameters.
#' @param log,log.p	logical; if TRUE, probabilities p are given as log(p).	
#' @param lower.tail logical; if TRUE (default), probabilities are 
#' P[X <= x], otherwise, P[X > x].
#' 
#' @seealso \link{IW}
#' 
#' @details 
#' The inverse weibull distribution with parameters \code{mu} and
#' \code{sigma} has density given by
#' 
#' \eqn{f(x) = \mu \sigma x^{-\sigma-1} \exp(\mu x^{-\sigma})}
#' 
#' for \eqn{x > 0}, \eqn{\mu > 0} and \eqn{\sigma > 0} 
#' 
#' @return 
#' \code{dIW} gives the density, \code{pIW} gives the distribution 
#' function, \code{qIW} gives the quantile function, \code{rIW}
#' generates random deviates and \code{hIW} gives the hazard function.
#'
#' @example examples/examples_dIW.R  
#'
#' @references
#' Almalki, S. J., & Nadarajah, S. (2014). Modifications of the 
#' Weibull distribution: A review. Reliability Engineering & 
#' System Safety, 124, 32-55.
#'
#' Drapella, A. (1993). The complementary Weibull distribution: 
#' unknown or just forgotten?. Quality and reliability engineering 
#' international, 9(4), 383-385.
#'
#' @export
dIW <- function(x, mu, sigma, log=FALSE){
  if (any(x < 0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(mu <= 0)) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  
  loglik <- log(mu*sigma) - (sigma+1)*log(x) - mu*(x^-sigma)
  
  if (log == FALSE) 
    density <- exp(loglik) 
  else density <- loglik
  return(density)  
}
#' @export
#' @rdname dIW
pIW <- function(q, mu, sigma, lower.tail=TRUE, log.p=FALSE){
  if (any(q < 0)) 
    stop(paste("q must be positive", "\n", ""))
  if (any(mu <= 0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
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
#' @rdname dIW
qIW <- function(p, mu, sigma, lower.tail = TRUE, log.p = FALSE){
  if (any(mu <= 0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  
  if (log.p == TRUE) 
    p <- exp(p)
  else p <- p
  if (lower.tail == TRUE) 
    p <- p
  else  p <- 1 - p
  if (any(p < 0) | any(p > 1)) 
    stop(paste("p must be between 0 and 1", "\n", ""))
  
  q <- ((-1/mu)*log(p))^(-1/sigma)
  q
}
#' @importFrom stats runif
#' @export
#' @rdname dIW
rIW <- function(n,mu,sigma){
  if(any(n <= 0))
    stop(paste("n must be positive","\n",""))
  if (any(mu <= 0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  
  n <- ceiling(n)
  p <- runif(n)
  r <- qIW(p, mu,sigma)
  r
}
#' @export
#' @rdname dIW
hIW<-function(x, mu, sigma){
  if (any(x < 0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(mu <= 0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  
  h <- dIW(x, mu, sigma, log=FALSE) / 
    pIW(q=x, mu, sigma, lower.tail=FALSE, log.p=FALSE)
  h
}

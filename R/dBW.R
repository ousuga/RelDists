#' The Beta Weibull distribution
#' 
#' @author Johan David Marin Benjumea, \email{johand.marin@@udea.edu.co}
#' 
#' @description 
#' Density, distribution function, quantile function, 
#' random generation and hazard function for the Beta Weibull distribution 
#' with parameters \code{mu}, \code{sigma}, \code{nu} and \code{tau}.
#' 
#' @param x,q	vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. 
#' @param mu parameter.
#' @param sigma parameter.
#' @param nu parameter.
#' @param tau parameter.
#' @param log,log.p	logical; if TRUE, probabilities p are given as log(p).	
#' @param lower.tail logical; if TRUE (default), probabilities are P[X <= x], otherwise, P[X > x].
#' 
#' @details 
#' The Beta Weibull Distribution with parameters \code{mu}, 
#' \code{sigma}, \code{nu} and \code{tau} has density given by
#' 
#' \eqn{f(x)= \frac{1}{B(\nu, \tau)} \mu \sigma x^{\sigma - 1} (1 - \exp(-\mu x^\sigma))^{\nu - 1} \exp(-\mu \tau x^\sigma),}
#' 
#' for \eqn{x > 0}, \eqn{\mu > 0}, \eqn{\sigma > 0}, \eqn{\nu > 0} and \eqn{\tau > 0}. 
#' 
#' @return 
#' \code{dBW} gives the density, \code{pBW} gives the distribution 
#' function, \code{qBW} gives the quantile function, \code{rBW}
#' generates random deviates and \code{hBW} gives the hazard function.
#'
#' @example examples/examples_dBW.R 
#' 
#' @references
#'\insertRef{almalki2014modifications}{RelDists}
#'
#'\insertRef{famoye2005beta}{RelDists}
#' 
#' @export
dBW <- function(x, mu, sigma, nu, tau, log=FALSE){
  if (any(x < 0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(mu <= 0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu <= 0)) 
    stop(paste("nu must be positive", "\n", ""))
  if (any(tau <= 0)) 
    stop(paste("tau must be positive", "\n", ""))
  
  
  loglik <- log(mu) + log(sigma) + (sigma - 1)*log(x) + 
    (nu - 1)*log(1 - exp(-mu*x^sigma)) - mu*tau*x^sigma - base::lbeta(nu, tau)
  
  if (log == FALSE) 
    density <- exp(loglik)
  else 
    density <- loglik
  return(density) 
}
#' @export
#' @rdname dBW
pBW <- function(q, mu, sigma, nu, tau, lower.tail=TRUE, log.p=FALSE){
  if (any(q < 0)) 
    stop(paste("q must be positive", "\n", ""))
  if (any(mu <= 0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu <= 0)) 
    stop(paste("nu must be positive", "\n", ""))
  if (any(tau <= 0)) 
    stop(paste("tau must be positive", "\n", ""))
  
  cdf  <- stats::pbeta(1 - exp(-mu*q^sigma), nu, tau)
  
  if (lower.tail == TRUE) cdf <- cdf
  else cdf <- 1 - cdf 
  if (log.p == FALSE) cdf <- cdf
  else cdf <- log(cdf)
  cdf
}
#' @export
#' @rdname dBW
qBW <- function(p, mu, sigma, nu, tau, lower.tail=TRUE, log.p=FALSE){
  if (any(mu <= 0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu <= 0)) 
    stop(paste("nu must be positive", "\n", ""))
  if (any(tau <= 0)) 
    stop(paste("tau must be positive", "\n", ""))
  
  if (log.p == TRUE) p <- exp(p)
  else p <- p
  if (lower.tail == TRUE) p <- p
  else p <- 1 - p
  if (any(p < 0) | any(p > 1)) 
    stop(paste("p must be between 0 and 1", "\n", ""))
  
  F.inv <- function(y, mu, sigma, nu, tau) {
    uniroot(function(x) {pBW(x, mu, sigma, nu, tau) - y},
            interval=c(0, 99999))$root
  }
  F.inv <- Vectorize(F.inv)
  F.inv(p, mu, sigma, nu, tau)
}
#' @importFrom stats runif
#' @export
#' @rdname dBW
rBW <- function(n, mu, sigma, nu, tau){
  if(any(n <= 0))
    stop(paste("n must be positive","\n",""))
  if (any(mu <= 0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu <= 0)) 
    stop(paste("nu must be positive", "\n", ""))
  if (any(tau <= 0)) 
    stop(paste("tau must be positive", "\n", ""))
  
  n <- ceiling(n)
  p <- runif(n)
  r <- qBW(p, mu, sigma, nu, tau)
  r
}
#' @export
#' @rdname dBW
hBW <- function(x, mu, sigma, nu, tau){
  if (any(x < 0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(mu <= 0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu <= 0)) 
    stop(paste("nu must be positive", "\n", ""))
  if (any(tau <= 0)) 
    stop(paste("tau must be positive", "\n", ""))
  
  h <- dBW(x, mu, sigma, nu, tau, log=FALSE) / 
    pBW(x, mu, sigma, nu, tau, lower.tail=FALSE, log.p=FALSE)
  h
}

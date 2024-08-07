#' The Kumaraswamy Inverse Weibull distribution
#' 
#' @author Freddy Hernandez, \email{fhernanb@@unal.edu.co}
#' 
#' @description 
#' Density, distribution function, quantile function, 
#' random generation and hazard function for the Kumaraswamy Inverse Weibull distribution 
#' with parameters \code{mu}, \code{sigma} and \code{nu}.
#' 
#' @param x,q	vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. 
#' @param mu parameter.
#' @param sigma parameter.
#' @param nu parameter.
#' @param log,log.p	logical; if TRUE, probabilities p are given as log(p).	
#' @param lower.tail logical; if TRUE (default), probabilities are P[X <= x], otherwise, P[X > x].
#' 
#' @details 
#' The Kumaraswamy Inverse Weibull Distribution with parameters \code{mu}, 
#' \code{sigma} and \code{nu} has density given by
#' 
#' \eqn{f(x)= \mu \sigma \nu x^{-\sigma - 1} \exp{- \mu x^{-\sigma}} (1 - \exp{- \mu x^{-\sigma}})^{\nu - 1},}
#' 
#' for \eqn{x > 0}, \eqn{\mu > 0}, \eqn{\sigma > 0} and \eqn{\nu > 0}. 
#' 
#' The KumIW distribution with \eqn{\nu=1} corresponds with the IW distribution.
#' 
#' @return 
#' \code{dKumIW} gives the density, \code{pKumIW} gives the distribution 
#' function, \code{qKumIW} gives the quantile function, \code{rKumIW}
#' generates random deviates and \code{hKumIW} gives the hazard function.
#'
#' @example  examples/examples_dKumIW.R
#' 
#' @references
#'\insertRef{almalki2014modifications}{RelDists}
#'
#'\insertRef{shahbaz2012kumaraswamy}{RelDists}
#' 
#' @export
dKumIW <- function(x, mu, sigma, nu, log=FALSE){
  if (any(x < 0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(mu <= 0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu <= 0 )) 
    stop(paste("nu must be positive", "\n", ""))
  
  loglik <- log(mu) + log(sigma) + log(nu) - (sigma+1)*log(x) - mu*x^(-sigma) +
    (nu-1)*log(1 - exp(-mu*x^(-sigma))) 
  
  if (log == FALSE) 
    density <- exp(loglik)
  else 
    density <- loglik
  return(density) 
}
#' @export
#' @rdname dKumIW
pKumIW <- function(q, mu, sigma, nu, lower.tail=TRUE, log.p=FALSE){
  if (any(q < 0)) 
    stop(paste("q must be positive", "\n", ""))
  if (any(mu <= 0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu <= 0)) 
    stop(paste("nu must be positive", "\n", ""))
  
  cdf  <- 1 - (1 - exp(-mu*q^(-sigma)))^nu
  
  if (lower.tail == TRUE) cdf <- cdf
  else cdf <- 1 - cdf 
  if (log.p == FALSE) cdf <- cdf
  else cdf <- log(cdf)
  cdf
}
#' @export
#' @rdname dKumIW
qKumIW <- function(p, mu, sigma, nu, lower.tail=TRUE, log.p=FALSE){
  if (any(mu <= 0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu <= 0)) 
    stop(paste("nu must be positive", "\n", ""))
  
  if (log.p == TRUE) p <- exp(p)
  else p <- p
  if (lower.tail == TRUE) p <- p
  else p <- 1 - p
  if (any(p < 0) | any(p > 1)) 
    stop(paste("p must be between 0 and 1", "\n", ""))
  q <- ((-1/mu) * log(1 - (1 - p)^(1/nu)))^(-1/sigma)
  q
}
#' @importFrom stats runif
#' @export
#' @rdname dKumIW
rKumIW <- function(n, mu, sigma, nu){
  if(any(n <= 0))
    stop(paste("n must be positive","\n",""))
  if (any(mu <= 0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu <= 0)) 
    stop(paste("nu must be positive", "\n", ""))
  
  n <- ceiling(n)
  p <- runif(n)
  r <- qKumIW(p, mu, sigma, nu)
  r
}
#' @export
#' @rdname dKumIW
hKumIW <- function(x, mu, sigma, nu){
  if (any(x < 0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(mu <= 0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu <= 0)) 
    stop(paste("nu must be positive", "\n", ""))
  
  h <- dKumIW(x, mu, sigma, nu, log=FALSE) / 
    pKumIW(x, mu, sigma, nu, lower.tail=FALSE, log.p=FALSE)
  h
}

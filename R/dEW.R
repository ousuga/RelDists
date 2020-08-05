#' The Exponentiated Weibull distribution
#' 
#' @description 
#' Density, distribution function, quantile function, 
#' random generation and hazard function for the exponentiated Weibull distribution with
#' parameters \code{mu}, \code{sigma} and \code{nu}.
#' 
#' @param x,q	vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. 
#' @param mu scale parameter.
#' @param sigma,nu shape parameters.
#' @param log,log.p	logical; if TRUE, probabilities p are given as log(p).	
#' @param lower.tail logical; if TRUE (default), probabilities are P[X <= x], otherwise, P[X > x].
#' 
#' @seealso \link{EW}
#' 
#' @details 
#' The Exponentiated Weibull Distribution with parameters \code{mu}, 
#' \code{sigma} and \code{nu} has density given by
#' 
#' \eqn{f(x)=\nu \mu \sigma x^{\sigma-1} \exp(-\mu x^\sigma) (1-\exp(-\mu x^\sigma))^{\nu-1},}
#' 
#' for \eqn{x > 0}, \eqn{\mu > 0}, \eqn{\sigma > 0} and \eqn{\nu > 0}. 
#' 
#' @return 
#' \code{dEW} gives the density, \code{pEW} gives the distribution 
#' function, \code{qEW} gives the quantile function, \code{rEW}
#' generates random deviates and \code{hEW} gives the hazard function.
#'
#' @example examples/examples_dEW.R
#' 
#' @export
dEW <- function(x, mu, sigma, nu, log = FALSE) {
  if (any(x < 0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(mu <= 0)) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu <= 0)) 
    stop(paste("nu must be positive", "\n", ""))
  loglik <- log(nu) + log(mu) + log(sigma) + (sigma - 1) * log(x) - mu * 
    (x^sigma) + (nu - 1) * log(1 - exp(-mu * (x^sigma)))
  if (log == FALSE) 
    density <- exp(loglik) 
  else density <- loglik
  return(density)
}
#' @export
#' @rdname dEW
pEW <- function(q, mu, sigma, nu, lower.tail = TRUE, log.p = FALSE) {
  if (any(q < 0)) 
    stop(paste("q must be positive", "\n", ""))
  if (any(mu <= 0)) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu <= 0)) 
    stop(paste("nu must be positive", "\n", ""))
  
  cdf <- (1 - exp(-mu * (q^sigma)))^nu
  
  if (lower.tail == TRUE) 
    cdf <- cdf 
  else cdf <- 1 - cdf
  if (log.p == FALSE) 
    cdf <- cdf 
  else cdf <- log(cdf)
  cdf
}
#' @export
#' @rdname dEW
qEW <- function(p, mu, sigma, nu, lower.tail = TRUE, log.p = FALSE) {
  if (any(mu <= 0)) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu <= 0)) 
    stop(paste("nu must be positive", "\n", ""))
  
  if (log.p == TRUE) 
    p <- exp(p)
  else p <- p
  if (lower.tail == TRUE) 
    p <- p 
  else p <- 1 - p
  if (any(p < 0) | any(p > 1)) 
    stop(paste("p must be between 0 and 1", "\n", ""))
  
  q <- ((-1/mu) * log(1 - p^(1/nu)))^(1/sigma)
  q
}
#' @importFrom stats runif
#' @export
#' @rdname dEW
rEW <- function(n, mu, sigma, nu) {
  if (any(n <= 0)) 
    stop(paste("n must be positive", "\n", ""))
  if (any(mu <= 0)) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu <= 0)) 
    stop(paste("nu must be positive", "\n", ""))
  
  n <- ceiling(n)
  p <- runif(n)
  r <- qEW(p, mu, sigma, nu)
  r
}
#' @export
#' @rdname dEW
hEW <- function(x, mu, sigma, nu) {
  if (any(x < 0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(mu <= 0)) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu <= 0)) 
    stop(paste("nu must be positive", "\n", ""))
  
  h <- dEW(x, mu, sigma, nu, log = FALSE) / 
    pEW(x, mu, sigma, nu, lower.tail=FALSE, log.p=FALSE)
  h
}
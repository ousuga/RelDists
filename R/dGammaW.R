#' The Gamma Weibull distribution
#' 
#' @author Johan David Marin Benjumea, \email{johand.marin@@udea.edu.co}
#' 
#' @description 
#' Density, distribution function, quantile function, 
#' random generation and hazard function for the Gamma Weibull distribution 
#' with parameters \code{mu}, \code{sigma}, \code{nu} and \code{tau}.
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
#' The Gamma Weibull Distribution with parameters \code{mu}, 
#' \code{sigma} and \code{nu} has density given by
#' 
#' \eqn{f(x)= \frac{\sigma \mu^{\nu}}{\Gamma(\nu)} x^{\nu \sigma - 1} \exp(-\mu x^\sigma),}
#' 
#' for \eqn{x > 0}, \eqn{\mu > 0}, \eqn{\sigma \geq 0} and \eqn{\nu > 0}. 
#' 
#' @return 
#' \code{dGammaW} gives the density, \code{pGammaW} gives the distribution 
#' function, \code{qGammaW} gives the quantile function, \code{rGammaW}
#' generates random deviates and \code{hGammaW} gives the hazard function.
#'
#' @examples  
#' ## The probability density function 
#' curve(dGammaW(x, mu = 0.5, sigma = 2, nu=1), from = 0, to = 10, 
#'       col = "red", las = 1, ylab = "The probability density function")
#' 
#' ## The cumulative distribution and the Reliability function
#' par(mfrow = c(1, 2))
#' curve(pGammaW(x, mu = 0.5, sigma = 2, nu=1), from = 0, to = 6, 
#' ylim = c(0, 1), col = "red", las = 1, ylab = "The cumulative distribution function")
#' curve(pGammaW(x, mu = 0.5, sigma = 2, nu=1, lower.tail = FALSE), 
#' from = 0, to = 6, ylim = c(0, 1), col = "red", las = 1, ylab = "The Reliability function")
#' 
#' ## The quantile function
#' p <- seq(from = 0, to = 0.99999, length.out = 100)
#' plot(x = qGammaW(p = p, mu = 0.5, sigma = 2, nu=1), y = p, 
#' xlab = "Quantile", las = 1, ylab = "Probability")
#' curve(pGammaW(x, mu = 0.5, sigma = 2, nu=1), from = 0, add = TRUE, 
#' col = "red")
#' 
#' ## The random function
#' hist(rGammaW(1000, mu = 0.5, sigma = 2, nu=1), freq = FALSE, xlab = "x", 
#' ylim = c(0, 1), las = 1, main = "")
#' curve(dGammaW(x, mu = 0.5, sigma = 2, nu=1),  from = 0, add = TRUE, 
#' col = "red", ylim = c(0, 1))
#' 
#' ## The Hazard function(
#' par(mfrow=c(1,1))
#' curve(hGammaW(x, mu = 0.5, sigma = 2, nu=1), from = 0, to = 2, 
#' ylim = c(0, 1), col = "red", ylab = "The hazard function", las = 1)
#' 
#' @export
dGammaW <- function(x, mu, sigma, nu, log=FALSE){
  if (any(x < 0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(mu <= 0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  loglik <- log(sigma) + nu*log(mu) + (nu*sigma - 1)*log(x) - mu*x^sigma -
    base::lgamma(nu)
  
  if (log == FALSE) 
    density <- exp(loglik)
  else 
    density <- loglik
  return(density) 
}
#' @export
#' @rdname dGammaW
pGammaW <- function(q, mu, sigma, nu,lower.tail=TRUE, log.p=FALSE){
  if (any(q < 0)) 
    stop(paste("q must be positive", "\n", ""))
  if (any(mu <= 0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu <= 0)) 
    stop(paste("nu must be positive", "\n", ""))
  
  cdf <- stats::pgamma(nu, mu*q^sigma, lower.tail=FALSE)
  
  if (lower.tail == TRUE) cdf <- cdf
  else cdf <- 1 - cdf 
  if (log.p == FALSE) cdf <- cdf
  else cdf <- log(cdf)
  cdf
}
#' @export
#' @rdname dGammaW
qGammaW <- function(p, mu, sigma, nu,
                    lower.tail=TRUE, log.p=FALSE){
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
  
  F.inv <- function(y, mu, sigma, nu) {
    uniroot(function(x) {pGammaW(x,mu,sigma,nu) - y},
            interval=c(0, 99999))$root
  }
  F.inv <- Vectorize(F.inv)
  F.inv(p, mu, sigma, nu)
}
#' @importFrom stats runif
#' @export
#' @rdname dGammaW
rGammaW <- function(n, mu, sigma, nu){
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
  r <- qGammaW(p, mu, sigma, nu)
  r
}
#' @export
#' @rdname dGammaW
hGammaW <- function(x, mu, sigma, nu){
  if (any(x < 0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(mu <= 0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu <= 0)) 
    stop(paste("nu must be positive", "\n", ""))
  
  h <- dGammaW(x, mu, sigma, nu, log=FALSE) / 
    pGammaW(x, mu, sigma, nu, lower.tail=FALSE, log.p=FALSE)
  h
}
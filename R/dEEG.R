#' The Extended Exponential Geometric distribution
#' 
#' @author Johan David Marin Benjumea, \email{johand.marin@@udea.edu.co}
#' 
#' @description 
#' Density, distribution function, quantile function, 
#' random generation and hazard function for the Extended Exponential Geometric distribution 
#' with parameters \code{mu} and \code{sigma}.
#' 
#' @param x,q	vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. 
#' @param mu parameter.
#' @param sigma parameter.
#' @param log,log.p	logical; if TRUE, probabilities p are given as log(p).	
#' @param lower.tail logical; if TRUE (default), probabilities are P[X <= x], otherwise, P[X > x].
#' 
#' @details 
#' The Extended Exponential Geometric distribution with parameters \code{mu}, 
#' and \code{sigma}has density given by
#' 
#' \eqn{f(x)= \mu \sigma \exp(-\mu x)(1 - (1 - \sigma)\exp(-\mu x))^{-2},}
#' 
#' for \eqn{x > 0}, \eqn{\mu > 0} and \eqn{\sigma > 0}. 
#' 
#' @return 
#' \code{dEEG} gives the density, \code{pEEG} gives the distribution 
#' function, \code{qEEG} gives the quantile function, \code{rEEG}
#' generates random deviates and \code{hEEG} gives the hazard function.
#'
#' @examples 
#' ## The probability density function 
#' par(mfrow=c(1,1))
#' curve(dEEG(x, mu = 1, sigma =3), from = 0, to = 10, 
#'       col = "red", las = 1, ylab = "The probability density function")
#' 
#' ## The cumulative distribution and the Reliability function
#' par(mfrow = c(1, 2))
#' curve(pEEG(x, mu = 1, sigma =3), from = 0, to = 10, 
#'       ylim = c(0, 1), col = "red", las = 1, ylab = "The cumulative distribution function")
#' curve(pEEG(x, mu = 1, sigma =3, lower.tail = FALSE), 
#'       from = 0, to = 6, ylim = c(0, 1), col = "red", las = 1, ylab = "The Reliability function")
#' 
#' ## The quantile function
#' p <- seq(from = 0, to = 0.99999, length.out = 100)
#' plot(x = qEEG(p = p, mu = 1, sigma =0.5), y = p, 
#'      xlab = "Quantile", las = 1, ylab = "Probability")
#' curve(pEEG(x, mu = 1, sigma =0.5), from = 0, add = TRUE, 
#'       col = "red")
#' 
#' ## The random function
#' hist(rEEG(1000, mu = 1, sigma =1), freq = FALSE, xlab = "x", 
#'      ylim = c(0, 0.9), las = 1, main = "")
#' curve(dEEG(x, mu = 1, sigma =1),  from = 0, add = TRUE, 
#'       col = "red", ylim = c(0, 0.8))
#' 
#' ## The Hazard function
#' par(mfrow=c(1,1))
#' curve(hEEG(x, mu = 1, sigma =0.5), from = 0, to = 2, 
#'       col = "red", ylab = "The hazard function", las = 1)
#' 
#' @references
#' \insertRef{almalki2014modifications}{RelDists}
#'
#' \insertRef{adamidis2005extension}{RelDists}
#' 
#' @export
dEEG <- function(x, mu, sigma, log=FALSE){
  if (any(x < 0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(mu <= 0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  
  loglik <- log(mu) + log(sigma) - mu*x - 2*log(1 - (1 - sigma)*exp(-mu*x))
  
  if (log == FALSE) 
    density <- exp(loglik)
  else 
    density <- loglik
  return(density) 
}
#' @export
#' @rdname dEEG
pEEG <- function(q, mu, sigma, lower.tail=TRUE, log.p=FALSE){
  if (any(q < 0)) 
    stop(paste("q must be positive", "\n", ""))
  if (any(mu <= 0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  
  cdf  <- (1 - exp(-mu*q))*(1-(1-sigma)*exp(-mu*q))^(-1)
  
  if (lower.tail == TRUE) cdf <- cdf
  else cdf <- 1 - cdf 
  if (log.p == FALSE) cdf <- cdf
  else cdf <- log(cdf)
  cdf
}
#' @export
#' @rdname dEEG
qEEG <- function(p, mu, sigma, lower.tail=TRUE, log.p=FALSE){
  if (any(mu <= 0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  
  if (log.p == TRUE) p <- exp(p)
  else p <- p
  if (lower.tail == TRUE) p <- p
  else p <- 1 - p
  if (any(p < 0) | any(p > 1)) 
    stop(paste("p must be between 0 and 1", "\n", ""))
  
  q <- -1/mu * log((1-p)/(1-p+sigma))
  q
}
#' @importFrom stats runif
#' @export
#' @rdname dEEG
rEEG <- function(n, mu, sigma){
  if(any(n <= 0))
    stop(paste("n must be positive","\n",""))
  if (any(mu <= 0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  
  n <- ceiling(n)
  p <- runif(n)
  r <- qEEG(p, mu, sigma)
  r
}
#' @export
#' @rdname dEEG
hEEG <- function(x, mu, sigma){
  if (any(x < 0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(mu <= 0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  
  h <- dEEG(x, mu, sigma, log=FALSE) / 
    pEEG(x, mu, sigma, lower.tail=FALSE, log.p=FALSE)
  h
}

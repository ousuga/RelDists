#' The Kumaraswamy Inverse Weibull distribution
#' 
#' @author Johan David Marin Benjumea, \email{johand.marin@@udea.edu.co}
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
#' \eqn{f(x)= \mu \sigma \nu x^{-\mu - 1} \exp{- \sigma x^{-\mu}} (1 - \exp{- \sigma x^{-\mu}})^{\nu - 1},}
#' 
#' for \eqn{x > 0}, \eqn{\mu > 0}, \eqn{\sigma > 0} and \eqn{\nu > 0}. 
#' 
#' @return 
#' \code{dKumIW} gives the density, \code{pKumIW} gives the distribution 
#' function, \code{qKumIW} gives the quantile function, \code{rKumIW}
#' generates random deviates and \code{hKumIW} gives the hazard function.
#'
#' @examples  
#' 
#' ## The probability density function 
#' par(mfrow = c(1, 1))
#' curve(dKumIW(x, mu = 1.5, sigma=  1.5, nu = 1), from = 0, to = 8.5, 
#'       col = "red", las = 1, ylab = "The probability density function")
#' 
#' ## The cumulative distribution and the Reliability function
#' par(mfrow = c(1, 2))
#' curve(pKumIW(x, mu = 1.5, sigma=  1.5, nu = 1), from = 0, to = 8.5, 
#'       ylim = c(0, 1), col = "red", las = 1, ylab = "The cumulative distribution function")
#' curve(pKumIW(x, mu = 1.5, sigma=  1.5, nu = 1, lower.tail = FALSE), 
#'       from = 0, to = 6, ylim = c(0, 1), col = "red", las = 1, ylab = "The Reliability function")
#' 
#' ## The quantile function
#' p <- seq(from = 0, to = 0.99999, length.out = 100)
#' plot(x = qKumIW(p=p, mu = 1.5, sigma=  1.5, nu = 10), y = p, 
#'      xlab = "Quantile", las = 1, ylab = "Probability")
#' curve(pKumIW(x, mu = 1.5, sigma=  1.5, nu = 10), from = 0, add = TRUE, 
#'       col = "red")
#' 
#' ## The random function
#' hist(rKumIW(1000, mu = 1.5, sigma=  1.5, nu = 5), freq = FALSE, xlab = "x", 
#'      las = 1, ylim = c(0, 1.5), main = "")
#' curve(dKumIW(x, mu = 1.5, sigma=  1.5, nu = 5), from = 0, to =8, add = TRUE, 
#'       col = "red")
#' 
#' ## The Hazard function
#' par(mfrow=c(1,1))
#' curve(hKumIW(x, mu = 1.5, sigma=  1.5, nu = 1), from = 0, to = 3, 
#'       ylim = c(0, 0.7), col = "red", ylab = "The hazard function", las = 1)
#' 
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
  
  
  loglik <- log(mu) + log(sigma) + log(nu) - (mu + 1)*log(x) - sigma*x^(-mu) +
    (nu -1)*log(1 - exp(- sigma*x^(-mu))) 
  
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
  
  cdf  <- 1 - (1 - exp(- sigma*q^(-mu)))^nu
  
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
  
  q <- (-sigma / (log(1 - (1 - p)^(1/nu))))
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

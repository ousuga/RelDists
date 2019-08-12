#' The Generalized Weibull Distribution
#' 
#' @author Johan David Marin Benjumea, \email{johand.marin@@udea.edu.co}
#' 
#' @description 
#' Density, distribution function, quantile function, 
#' random generation and hazard function for the Generalized Weibull distribution 
#' with parameters \code{mu}, \code{sigma}, and \code{nu}.
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
#' The Generalized Weibull Distribution with parameters \code{mu}, 
#' \code{sigma} and \code{nu} has density given by
#' 
#' \eqn{f(x)= \mu \sigma x^{\sigma - 1} (1 - \mu \nu x^\sigma)^{\frac{1}{\nu} - 1},}
#' 
#' for \eqn{x > 0}, \eqn{\mu > 0}, \eqn{\sigma > 0} and \eqn{0 < \nu < 1}. 
#' 
#' @return 
#' \code{dGWF} gives the density, \code{pGWF} gives the distribution 
#' function, \code{qGWF} gives the quantile function, \code{rGWF}
#' generates random deviates and \code{hGWF} gives the hazard function.
#'
#' @examples 
#' 
#' ## The probability density function 
#' par(mfrow = c(1, 1))
#' curve(dGWF(x, mu = 0.01, sigma = 2.9, nu=0.2), from = 0, to = 8.5, 
#'       col = "red", las = 1, ylab = "The probability density function")
#' 
#' ## The cumulative distribution and the Reliability function
#' par(mfrow = c(1, 2))
#' curve(pGWF(x, mu = 0.01, sigma = 2.9, nu=0.2), from = 0, to = 8.5, 
#'       ylim = c(0, 1), col = "red", las = 1, ylab = "The cumulative distribution function")
#' curve(pGWF(x, mu = 0.01, sigma = 2.9, nu=0.2, lower.tail = FALSE), 
#'       from = 0, to = 6, ylim = c(0, 1), col = "red", las = 1, ylab = "The Reliability function")
#' 
#' ## The quantile function
#' p <- seq(from = 0, to = 0.99999, length.out = 100)
#' plot(x = qGWF(p=p, mu = 0.01, sigma = 2.9, nu=0.2), y = p, 
#'      xlab = "Quantile", las = 1, ylab = "Probability")
#' curve(pGWF(x, mu = 0.01, sigma = 2.9, nu=0.2), from = 0, add = TRUE, 
#'       col = "red")
#' 
#' ## The random function
#' hist(rGWF(1000, mu = 0.01, sigma = 2.9, nu=0.2), freq = FALSE, xlab = "x", 
#'      las = 1, ylim = c(0, 0.3), main = "")
#' curve(dGWF(x, mu = 0.01, sigma = 2.9, nu=0.2), from = 0, to =8, add = TRUE, 
#'       col = "red")
#' 
#' ## The Hazard function
#' par(mfrow=c(1,1))
#' curve(hGWF(x, mu = 0.01, sigma = 2.9, nu=0.2), from = 0, to = 3, 
#'       ylim = c(0, 0.26), col = "red", ylab = "The hazard function", las = 1)
#' 
#' @references
#'\insertRef{almalki2014modifications}{RelDists}
#'
#'\insertRef{mudholkar1994generalized}{RelDists}
#' 
#' @export
dGWF <- function(x, mu, sigma, nu, log=FALSE){
  if (any(x < 0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(mu <= 0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu <= 0  | nu >= 1  )) 
    stop(paste("nu must be between zero and one", "\n", ""))
  
  loglik <- log(mu) + log(sigma) + (sigma - 1)*log(x) + 
    (1/nu - 1)*log(1 - mu*nu*x^sigma)
  
  if (log == FALSE) 
    density <- exp(loglik)
  else 
    density <- loglik
  return(density) 
}
#' @export
#' @rdname dGWF
pGWF <- function(q, mu, sigma, nu, lower.tail=TRUE, log.p=FALSE){
  if (any(q < 0)) 
    stop(paste("q must be positive", "\n", ""))
  if (any(mu <= 0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu <= 0  | nu >= 1  )) 
    stop(paste("nu must be between zero and one", "\n", ""))
  
  cdf  <- 1 - (1 - mu*nu*q^sigma)^(1/nu)
  
  if (lower.tail == TRUE) cdf <- cdf
  else cdf <- 1 - cdf 
  if (log.p == FALSE) cdf <- cdf
  else cdf <- log(cdf)
  cdf
}
#' @export
#' @rdname dGWF
qGWF <- function(p, mu, sigma, nu, lower.tail=TRUE, log.p=FALSE){
  if (any(mu <= 0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu <= 0  | nu >= 1  )) 
    stop(paste("nu must be between zero and one", "\n", ""))
  
  if (log.p == TRUE) p <- exp(p)
  else p <- p
  if (lower.tail == TRUE) p <- p
  else p <- 1 - p
  if (any(p < 0) | any(p > 1)) 
    stop(paste("p must be between 0 and 1", "\n", ""))
  
  q <- ((1-(1-p)^nu)/(mu*nu))^(1/sigma)
  q
}
#' @importFrom stats runif
#' @export
#' @rdname dGWF
rGWF <- function(n, mu, sigma, nu){
  if(any(n <= 0))
    stop(paste("n must be positive","\n",""))
  if (any(mu <= 0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu <= 0  | nu >= 1  )) 
    stop(paste("nu must be between zero and one", "\n", ""))
  
  n <- ceiling(n)
  p <- runif(n)
  r <- qGWF(p, mu, sigma, nu)
  r
}
#' @export
#' @rdname dGWF
hGWF <- function(x, mu, sigma, nu){
  if (any(x < 0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(mu <= 0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu <= 0  | nu >= 1  )) 
    stop(paste("nu must be between zero and one", "\n", ""))
  
  h <- dGWF(x, mu, sigma, nu, log=FALSE) / 
    pGWF(x, mu, sigma, nu, lower.tail=FALSE, log.p=FALSE)
  h
}
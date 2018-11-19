#' The Exponentiated Weibull Distribution
#' 
#' @description 
#' Density, distribution function, quantile function, 
#' random generation  and hazard function for the exponentiated weibull distribution with
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
#' #' @details 
#' The Exponentiated Weibull Distribution with parameters \code{mu}, 
#' \code{sigma} and \code{nu} has density given by
#' 
#' \eqn{f(x)=\nu \mu \sigma x^{\sigma-1} \exp(-\mu x^\sigma) (1-\exp(-\mu x^\sigma))^{\nu-1},}
#' 
#' for x > 0. 
#' 
#' @return 
#' \code{dEW} gives the density, \code{pEW} gives the distribution 
#' function, \code{qEW} gives the quantile function, \code{rEW}
#' generates random deviates and \code{hEW} gives the hazard function.
#'
#' @examples  
#' ## The probability density function
#' curve(dEW(x, mu=2, sigma=1.5, nu=0.5), from=0, to=2,
#'       ylim=c(0, 2.5), col="red", las=1, ylab="The probability density function") 
#' 
#' ## The cumulative distribution and the Reliability function
#' par(mfrow=c(1, 2))
#' curve(pEW(x, mu=2, sigma=1.5, nu=0.5), from=0, to=2,  col="red",
#'       las=1, ylab="The cumulative distribution function")
#' curve(pEW(x, mu=2, sigma=1.5, nu=0.5, lower.tail=FALSE), from=0,
#'       to=2, col="red", las=1, ylab="The Reliability function")
#' 
#' ## The quantile function
#' p <- seq(from=0, to=0.99999, length.out=100)
#' plot(x=qEW(p, mu=2, sigma=1.5, nu=0.5), y=p, xlab="Quantile", las=1, ylab="Probability")
#' curve(pEW(x, mu=2, sigma=1.5, nu=0.5),  from=0, add=TRUE, col="red")
#' 
#' ## The random function
#' hist(rEW(n=10000, mu=2, sigma=1.5, nu=0.5), freq=FALSE, xlab="x", las=1, main="")
#' curve(dEW(x, mu=2, sigma=1.5, nu=0.5),  from=0, add=TRUE, col="red") 
#' 
#' ## The Hazard function
#' curve(hEW(x, mu=2,sigma=1.5, nu=0.5), from=0, to=2, ylim=c(0, 7), col="red",
#'       ylab="The Hazard function")
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
    density <- exp(loglik) else density <- loglik
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
    cdf <- cdf else cdf <- 1 - cdf
  if (log.p == FALSE) 
    cdf <- cdf else cdf <- log(cdf)
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
    p <- exp(p) else p <- p
    if (lower.tail == TRUE) 
      p <- p else p <- 1 - p
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
    pEW(q = x, mu, sigma, nu, lower.tail = FALSE, log.p = FALSE)
  h
}


#' The Kumaraswamy Weibull Distribution
#' 
#' @description 
#' Density, distribution function, quantile function, 
#' random generation and hazard function for the Kumaraswamy Weibull distribution 
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
#' The Kumaraswamy Weibull Distribution with parameters \code{mu}, 
#' \code{sigma}, \code{nu} and \code{tau} has density given by
#' 
#' \eqn{f(x)=\nu \tau \mu \sigma x^{\sigma - 1} \exp^{-\mu x ^\sigma} [1 - \exp^{-\mu x ^ \sigma}]^{\nu-1}[1-(1-\exp^{-\mu x ^ \sigma})^\nu]^{\tau-1},}
#' 
#' for x > 0. 
#' 
#' @return 
#' \code{dKW} gives the density, \code{pKW} gives the distribution 
#' function, \code{qKW} gives the quantile function, \code{rKW}
#' generates random deviates and \code{hKW} gives the hazard function.
#'
#' @examples  
#' ## The probability density function
#' curve(dKW(x, mu=3, sigma=0.8, nu=2.0, tau=1.5), from=0, to=2,
#'       ylim=c(0, 2.5), col="red", las=1, ylab="f(x)")
#' 
#' ## The cumulative distribution and the Reliability function
#' par(mfrow=c(1, 2))
#' curve(pKW(x, mu=3, sigma=0.8, nu=2.0, tau=1.5),
#'       from=0, to=2,  col="red", las=1, ylab="F(x)")
#' curve(pKW(x, mu=3, sigma=0.8, nu=2.0, tau=1.5, lower.tail=FALSE),
#'       from=0, to=2, col="red", las=1, ylab="S(x)")
#' 
#' ## The quantile function
#' p <- seq(from=0, to=0.99999, length.out=100)
#' plot(x=qKW(p, mu=3, sigma=0.8, nu=2.0, tau=1.5), y=p, xlab="Quantile",
#'      las=1, ylab="Probability")
#' curve(pKW(x, mu=3, sigma=0.8, nu=2.0, tau=1.5), from=0, add=TRUE, col="red")
#' 
#' ## The random function
#' hist(rKW(n=10000, mu=3, sigma=0.8, nu=2.0, tau=1.5), freq=FALSE,
#'      xlab="x", las=1, main="")
#' curve(dKW(x, mu=3, sigma=0.8, nu=2.0, tau=1.5), from=0, add=TRUE, col="red")
#' 
#' ## The Hazard function
#' curve(hKW(x, mu=3, sigma=0.8, nu=2.0, tau=1.5), from=0, to=2, ylim=c(0, 7),
#'       col="red", ylab="Hazard function", las=1)
#'
#' @export
dKW <- function(x, mu, sigma, nu, tau, log=FALSE){
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
  term <- -mu * (x^sigma)
  loglik <- log(nu*tau*mu*sigma) + (sigma-1)*log(x) + term + 
    (nu-1)*log(1-exp(term)) + (tau-1)*log(1-(1-exp(term))^nu)
  if (log == FALSE) 
    density <- exp(loglik)
  else 
    density <- loglik
  return(density) 
}
#' @export
#' @rdname dKW
pKW <- function(q, mu, sigma, nu, tau, 
                lower.tail=TRUE, log.p=FALSE){
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
  cdf <- 1 - (1 - (1 - exp(-mu*q^sigma))^nu)^tau
  if (lower.tail == TRUE) cdf <- cdf
  else cdf <- 1 - cdf 
  if (log.p == FALSE) cdf <- cdf
  else cdf <- log(cdf)
  cdf
}
#' @export
#' @rdname dKW
qKW <- function(p, mu, sigma, nu, tau, 
                lower.tail=TRUE, log.p=FALSE){
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
  q <- ((-1/mu)*(log(1-(1-(1-p)^(1/tau))^(1/nu))))^(1/sigma)
  q
}
#' @importFrom stats runif
#' @export
#' @rdname dKW
rKW <- function(n, mu, sigma, nu, tau){
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
  r <- qKW(p, mu, sigma, nu, tau)
  r
}
#' @export
#' @rdname dKW
hKW <- function(x, mu, sigma, nu, tau){
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
  h <- dKW(x, mu, sigma, nu, tau, log=FALSE) / 
    pKW(x, mu, sigma, nu, tau, lower.tail=FALSE, log.p=FALSE)
  h
}

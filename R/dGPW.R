#' The Generalized power Weibull Distribution
#' 
#' @description 
#' Density, distribution function, quantile function, 
#' random generation and hazard function for the Generalized power Weibull Distribution with
#' parameters \code{mu}, \code{sigma} and \code{nu}.
#' 
#' @param x,q	vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. 
#' @param mu parameter.
#' @param sigma parameters.
#' @param nu parameters.
#' @param log,log.p	logical; if TRUE, probabilities p are given as log(p).	
#' @param lower.tail logical; if TRUE (default), probabilities are P[X <= x], otherwise, P[X > x].
#' 
#' @details 
#' The Generalized power Weibull Distribution with parameters \code{mu}, 
#' \code{sigma} and \code{nu} has density given by
#' 
#' \eqn{f(x)=\mu\sigma\nu^{-1} x^{\sigma-1}(1+\mu x^\sigma)^{\frac 1{\nu}-1}\;exp{[1-(1+\mu x^\sigma)^{\frac 1{\nu}}]}}
#' 
#' for x > 0. 
#' 
#' @return 
#' \code{dGPW} gives the density, \code{pGPW} gives the distribution 
#' function, \code{qGPW} gives the quantile function, \code{rGPW}
#' generates random deviates and \code{hGPW} gives the hazard function.
#'
#' @examples  
#'
#'## The probability density function
#'curve(dGPW(x, mu=0.2, sigma=2.0, nu=0.3), 
#'      xlim=c(0,2), ylim=c(0, 2.5), col="red", las=1, ylab="f(x)")
#'
#'
#'## The cumulative distribution and the Reliability function
#'par(mfrow=c(1, 2))
#'curve(pGPW(x, mu=0.2, sigma=2.0, nu=0.3), 
#'      xlim=c(0,2), col="red", las=1, ylab="F(x)")
#'curve(pGPW(x, mu=2, sigma=1.5, nu=0.5, lower.tail=FALSE),
#'      from=0, to=2, col="red", las=1, ylab="R(x)")
#'
#'## The quantile function
#'p <- seq(from=0, to=0.99999, length.out=100)
#'plot(x=qGPW(p, mu=0.2, sigma=2.0, nu=0.3), y=p, xlab="Quantile",
#'     las=1, ylab="Probability")
#'curve(pGPW(x, mu=0.2, sigma=2.0, nu=0.3), xlim=c(0,2), add=TRUE, col="red")
#'
#'## The random function
#'hist(rGPW(n=10000, mu=0.2, sigma=2.0, nu=0.3), freq=FALSE,
#'     xlab="x", las=1, main="")
#'curve(dGPW(x, mu=0.2, sigma=2.0, nu=0.3), xlim=c(0,2), add=TRUE, col="red")
#'
#'## The Hazard function
#'curve(hGPW(x, mu=0.2, sigma=2.0, nu=0.3), xlim=c(0,2), ylim=c(0, 7),
#'      col="red", ylab="Hazard function", las=1)
#'
#'
#' @export
dGPW <- function(x, mu, sigma, nu, log=FALSE){
  if (any(x < 0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(mu <= 0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu <= 0)) 
    stop(paste("nu must be positive", "\n", ""))
  loglik <- log(mu) + log(sigma) - log(nu) + (sigma-1) * log(x) +
    ((1/nu)-1)*log(1+mu*(x^sigma)) + (1-(1+mu*(x^sigma))^(1/nu))
  if (log == FALSE) 
    density <- exp(loglik)
  else 
    density <- loglik
  return(density)
}
#' @export
#' @rdname dGPW
pGPW <- function(q, mu, sigma, nu, lower.tail=TRUE, log.p = FALSE){
  if (any(q < 0)) 
    stop(paste("q must be positive", "\n", ""))
  if (any(mu <= 0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu <= 0)) 
    stop(paste("nu must be positive", "\n", ""))
  cdf <- 1 - exp(1 - (1 + mu*(q^sigma))^(1/nu))
  if (lower.tail == TRUE) 
    cdf <- cdf
  else cdf <- 1 - cdf
  if (log.p == FALSE) 
    cdf <- cdf
  else cdf <- log(cdf)
  cdf 
}
#' @export
#' @rdname dGPW
qGPW <- function(p, mu, sigma, nu, lower.tail=TRUE, log.p=FALSE){
  if (any(mu <= 0 )) 
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
  else p <- 1-p
  if (any(p < 0) | any(p > 1)) 
    stop(paste("p must be between 0 and 1", "\n", ""))
  term <- 1-log(1-p)
  q <- (((term^nu)-1)/mu)^(1/sigma)
  q
}
#' @importFrom stats runif
#' @export
#' @rdname dGPW
rGPW <- function(n, mu, sigma, nu){
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
  r <- qGPW(p, mu, sigma, nu)
  r
}
#' @export
#' @rdname dGPW
hGPW<-function(x, mu, sigma, nu){
  if (any(x<0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(mu <= 0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu <= 0)) 
    stop(paste("nu must be positive", "\n", ""))
  h <- dGPW(x, mu, sigma, nu, log=FALSE)/pGPW(q=x, mu, sigma, nu, lower.tail=FALSE, log.p=FALSE)
  h
}


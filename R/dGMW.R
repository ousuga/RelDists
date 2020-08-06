#' The Generalized modified Weibull distribution 
#' 
#' @description 
#' Density, distribution function, quantile function, 
#' random generation and hazard function for the generalized 
#' modified weibull distribution with parameters \code{mu}, 
#' \code{sigma}, \code{nu} and \code{tau}.
#' 
#' @param x,q	vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. 
#' @param mu scale parameter.
#' @param sigma shape parameter.
#' @param nu shape parameter.
#' @param tau acceleration parameter.
#' @param log,log.p	logical; if TRUE, probabilities p are given as log(p).	
#' @param lower.tail logical; if TRUE (default), probabilities are 
#' P[X <= x], otherwise, P[X > x].
#' 
#' @details
#' The generalized modified weibull with parameters \code{mu}, 
#' \code{sigma}, \code{nu} and \code{tau} has density given by
#' 
#' \eqn{f(x)= \mu \sigma x^{\nu - 1}(\nu + \tau x) \exp(\tau x - \mu x^{\nu} e^{\tau x})
#' [1 - \exp(- \mu x^{\nu} e^{\tau x})]^{\sigma-1},}
#' 
#' for x>0.
#' 
#' @return 
#' \code{dGMW} gives the density, \code{pGMW} gives the distribution 
#' function, \code{qGMW} gives the quantile function, \code{rGMW}
#' generates random deviates and \code{hGMW} gives the hazard function.
#' 
#' @example examples/examples_dGMW.R
#'
#' @export
dGMW <- function(x, mu, sigma, nu, tau, log=FALSE){
  if (any(x<0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(mu<=0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma<=0)) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu<0)) 
    stop(paste("nu must be positive", "\n", ""))
  if (any(tau<0)) 
    stop(paste("tau must be positive", "\n", ""))
  
  loglik <- log(mu*sigma) + (nu-1)*log(x) + log(nu + tau*x) +
            tau*x - mu*(x^nu)*exp(tau*x)  +
            (sigma-1)*log(1-exp(-mu*(x^nu)*exp(tau*x) ))
  
  if (log==FALSE) 
    density <- exp(loglik)
  else density <- loglik
  return(density)
}
#' @export
#' @rdname dGMW
pGMW <- function(q, mu, sigma, nu, tau, lower.tail=TRUE, log.p=FALSE){
  if (any(q<0)) 
    stop(paste("q must be positive", "\n", ""))
  if (any(mu<=0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma<=0)) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu<0)) 
    stop(paste("nu must be positive", "\n", ""))
  if (any(tau<0)) 
    stop(paste("tau must be positive", "\n", ""))
  
  cdf <- (1-exp(-mu*(q^nu)*exp(tau*q)))^sigma
  
  if (lower.tail==TRUE) 
    cdf <- cdf
  else cdf <- 1 - cdf
  if (log.p==FALSE) 
    cdf <- cdf
  else cdf <- log(cdf)
  cdf
}
#' @export
#' @rdname dGMW
qGMW <- function(p, mu, sigma, nu, tau, lower.tail=TRUE, log.p=FALSE) {
  if (any(mu<=0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma<=0)) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu<0)) 
    stop(paste("nu must be positive", "\n", ""))
  if (any(tau<0)) 
    stop(paste("tau must be positive", "\n", ""))
  
  if (log.p==TRUE) 
    p <- exp(p)
  else p <- p
  if (lower.tail==TRUE) 
    p <- p
  else p <- 1 - p
  if (any(p < 0) | any(p > 1)) 
    stop(paste("p must be between 0 and 1", "\n", ""))
  
  cdf <- function(x, mu, sigma, nu, tau){
    (1-exp((-mu*(x^nu))*exp(tau*x)))^sigma
  }
  cdf1 <- function(x, mu, sigma, nu, tau, p) {cdf(x, mu, sigma, nu, tau) - p}
    root.cdf1 <- function(mu, sigma, nu, tau, p) {
    uniroot(cdf1, interval=c(0,1e+06), mu, sigma, nu, tau, p)$root
  }
  root.cdf1 <- Vectorize(root.cdf1)
  q <- root.cdf1(mu, sigma, nu, tau, p)
  q
}
#' @importFrom stats runif
#' @export
#' @rdname dGMW
rGMW <- function(n, mu, sigma, nu, tau){
  if (any(mu<=0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma<=0)) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu<0)) 
    stop(paste("nu must be positive", "\n", ""))
  if (any(tau<0)) 
    stop(paste("tau must be positive", "\n", ""))
  
  n <- ceiling(n)
  p <- runif(n)
  r <- qGMW(p, mu, sigma, nu, tau)
  r
}
#' @export
#' @rdname dGMW
hGMW <- function(x, mu, sigma, nu, tau, log=FALSE){
  if (any(x<0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(mu<=0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma<=0)) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu<0)) 
    stop(paste("nu must be positive", "\n", ""))
  if (any(tau<0)) 
    stop(paste("tau must be positive", "\n", ""))
  
  h <- dGMW(x, mu, sigma, nu, tau, log=FALSE) / 
       pGMW(q=x, mu, sigma, nu, tau, lower.tail=FALSE, log.p=FALSE)
  h
}
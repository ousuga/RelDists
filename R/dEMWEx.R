#' The Exponentiated Modified Weibull Extension distribution
#' 
#' @author Johan David Marin Benjumea, \email{johand.marin@@udea.edu.co}
#' 
#' @description 
#' Density, distribution function, quantile function, 
#' random generation and hazard function for the Exponentiated Modifien Weibull Extension distribution 
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
#' @param lower.tail logical; if TRUE (default), probabilities are 
#' P[X <= x], otherwise, P[X > x].
#' 
#' @seealso \link{EMWEx}
#' 
#' @details 
#' The Exponentiated Modified Weibull Extension Distribution with parameters \code{mu}, 
#' \code{sigma}, \code{nu} and \code{tau} has density given by
#' 
#' \eqn{f(x)= \nu \sigma \tau (\frac{x}{\mu})^{\sigma-1} \exp((\frac{x}{\mu})^\sigma +
#' \nu \mu (1- \exp((\frac{x}{\mu})^\sigma))) 
#' (1 - \exp (\nu\mu (1- \exp((\frac{x}{\mu})^\sigma))))^{\tau-1} ,}
#' 
#' for \eqn{x > 0}, \eqn{\nu> 0}, \eqn{\mu > 0}, \eqn{\sigma> 0} and \eqn{\tau > 0}. 
#' 
#' @return 
#' \code{dEMWEx} gives the density, \code{pEMWEx} gives the distribution 
#' function, \code{qEMWEx} gives the quantile function, \code{rEMWEx}
#' generates random deviates and \code{hEMWEx} gives the hazard function.
#'
#' @example examples/examples_dEMWEx.R
#' 
#' @references
#' Almalki, S. J., & Nadarajah, S. (2014). Modifications of the 
#' Weibull distribution: A review. Reliability Engineering & 
#' System Safety, 124, 32-55.
#'
#' Sarhan, A. M., & Apaloo, J. (2013). Exponentiated modified 
#' Weibull extension distribution. Reliability Engineering & 
#' System Safety, 112, 137-144.
#' 
#' @export
dEMWEx <- function(x, mu, sigma, nu, tau, log=FALSE){
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
  
  exp1 <- (x/mu)^sigma
  exp2 <- mu*nu*(1-exp(exp1))
  lik <- nu*sigma*tau * (x/mu)^(sigma - 1) * exp(exp1 + exp2)* (1- exp(exp2))^(tau-1)
  #  loglik <- log(sigma) + log(nu) + log(tau) + (sigma - 1)*(log(x/mu))+ exp1 + exp2 + (tau - 1)*log(1- exp(exp2))
  
  if (log == FALSE) 
    density <- lik
  else 
    density <- log(lik)
  return(density) 
}
#' @export
#' @rdname dEMWEx
pEMWEx <- function(q, mu, sigma, nu, tau, lower.tail=TRUE, log.p=FALSE){
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
  
  exp1 <- (q/mu)^sigma
  exp2 <- mu*nu*(1-exp(exp1))
  cdf  <- (1- exp(exp2))^tau
  
  if (lower.tail == TRUE) cdf <- cdf
  else cdf <- 1 - cdf 
  if (log.p == FALSE) cdf <- cdf
  else cdf <- log(cdf)
  cdf
}
#' @export
#' @rdname dEMWEx
qEMWEx <- function(p, mu, sigma, nu, tau, lower.tail=TRUE, log.p=FALSE){
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
  
  q <- mu*(log(1-((log(1-p^(1/tau)))/(nu*mu))))^(1/sigma)
  q
}
#' @importFrom stats runif
#' @export
#' @rdname dEMWEx
rEMWEx <- function(n, mu, sigma, nu, tau){
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
  r <- qEMWEx(p, mu, sigma, nu, tau)
  r
}
#' @export
#' @rdname dEMWEx
hEMWEx <- function(x, mu, sigma, nu, tau){
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
  
  h <- dEMWEx(x, mu, sigma, nu, tau, log=FALSE) / 
    pEMWEx(x, mu, sigma, nu, tau, lower.tail=FALSE, log.p=FALSE)
  h
}

#' The Beta Generalized Exponentiated  distribution
#' 
#' @author Johan David Marin Benjumea, \email{johand.marin@@udea.edu.co}
#' 
#' @description 
#' Density, distribution function, quantile function, 
#' random generation and hazard function for the Beta Generalized Exponentiated  distribution 
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
#' @param lower.tail logical; if TRUE (default), probabilities 
#' are P[X <= x], otherwise, P[X > x].
#' 
#' @seealso \link{BGE}
#' 
#' @details 
#' The Beta Generalized Exponentiated  Distribution with parameters \code{mu}, 
#' \code{sigma}, \code{nu} and \code{tau} has density given by
#' 
#' \eqn{f(x)= \frac{\nu \tau}{B(\mu, \sigma)} \exp(-\nu x)(1- \exp(-\nu x))^{\tau \mu - 1} (1 - (1- \exp(-\nu x))^\tau)^{\sigma -1},}
#' 
#' for \eqn{x > 0}, \eqn{\mu > 0}, \eqn{\sigma > 0}, \eqn{\nu > 0} and \eqn{\tau > 0}. 
#' 
#' @return 
#' \code{dBGE} gives the density, \code{pBGE} gives the distribution 
#' function, \code{qBGE} gives the quantile function, \code{rBGE}
#' generates random deviates and \code{hBGE} gives the hazard function.
#'
#' @example examples/examples_dBGE.R    
#' 
#' @references
#' Almalki, S. J., & Nadarajah, S. (2014). Modifications of the 
#' Weibull distribution: A review. Reliability Engineering & 
#' System Safety, 124, 32-55.
#'
#' Barreto-Souza, W., Santos, A. H., & Cordeiro, G. M. (2010). 
#' The beta generalized exponential distribution. Journal of 
#' statistical Computation and Simulation, 80(2), 159-172.
#' 
#' @export
dBGE <- function(x, mu, sigma, nu, tau, log=FALSE){
  if (any(x < 0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(mu <= 0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu <= 0 )) 
    stop(paste("nu must be positive", "\n", ""))
  if (any(tau <= 0)) 
    stop(paste("tau must be positive", "\n", ""))
  
  exp1 <- 1 - exp(-nu*x)
  loglik <- log(tau) + log(nu) - nu*x + (tau*mu - 1)*log(exp1) + 
    (sigma - 1)*log(1- exp1^tau) - base::lbeta(mu, sigma)
  
  if (log == FALSE) 
    density <- exp(loglik)
  else 
    density <- loglik
  return(density) 
}
#' @export
#' @rdname dBGE
pBGE <- function(q, mu, sigma, nu, tau, lower.tail=TRUE, log.p=FALSE){
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
  
  exp1 <- (1 - exp(-nu*q))^tau
  cdf  <- stats::pbeta(exp1, mu, sigma)
  
  if (lower.tail == TRUE) cdf <- cdf
  else cdf <- 1 - cdf 
  if (log.p == FALSE) cdf <- cdf
  else cdf <- log(cdf)
  cdf
}
#' @export
#' @rdname dBGE
qBGE <- function(p, mu, sigma, nu, tau, lower.tail=TRUE, log.p=FALSE){
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
  
  F.inv <- function(y, mu, sigma, nu, tau) {
    uniroot(function(x) {pBGE(x, mu, sigma, nu, tau) - y},
            interval=c(0, 99999))$root
  }
  F.inv <- Vectorize(F.inv)
  F.inv(p, mu, sigma, nu, tau)
}
#' @importFrom stats runif
#' @export
#' @rdname dBGE
rBGE <- function(n, mu, sigma, nu, tau){
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
  r <- qBGE(p, mu, sigma, nu, tau)
  r
}
#' @export
#' @rdname dBGE
hBGE <- function(x, mu, sigma, nu, tau){
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
  
  h <- dBGE(x, mu, sigma, nu, tau, log=FALSE) / 
    pBGE(x, mu, sigma, nu, tau, lower.tail=FALSE, log.p=FALSE)
  h
}

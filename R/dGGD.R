#' The Generalized Gompertz  distribution
#' 
#' @author Johan David Marin Benjumea, \email{johand.marin@@udea.edu.co}
#' 
#' @description 
#' Density, distribution function, quantile function, 
#' random generation and hazard function for the generalized Gompertz distribution with
#' parameters \code{mu} \code{sigma} and \code{nu}.
#' 
#' @param x,q	vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. 
#' @param mu,nu scale parameter.
#' @param sigma shape parameters.
#' @param log,log.p	logical; if TRUE, probabilities p are given as log(p).	
#' @param lower.tail logical; if TRUE (default), probabilities are P[X <= x], otherwise, P[X > x].
#' 
#' @details 
#' The Generalized Gompertz  Distribution with parameters \code{mu}, 
#' \code{sigma} and \code{nu} has density given by
#' 
#' \eqn{f(x)= \nu \mu \exp(-\frac{\mu}{\sigma}(\exp(\sigma x - 1))) (1 - \exp(-\frac{\mu}{\sigma}(\exp(\sigma x - 1))))^{(\nu - 1)} ,}
#' 
#' for \eqn{x \geq 0}, \eqn{\mu > 0}, \eqn{\sigma \geq 0} and \eqn{\nu > 0}. 
#' 
#' @return 
#' \code{dGGD} gives the density, \code{pGGD} gives the distribution 
#' function, \code{qGGD} gives the quantile function, \code{rGGD}
#' generates random deviates and \code{hGGD} gives the hazard function.
#'
#' @example examples/examples_dGGD.R
#' 
#' @references
#' El-Gohary, A., Alshamrani, A., & Al-Otaibi, A. N. (2013). 
#' The generalized Gompertz distribution. Applied mathematical 
#' modelling, 37(1-2), 13-24.
#'
#' @export
dGGD <- function(x, mu, sigma, nu, log=FALSE){
  if (any(x < 0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(mu <= 0)) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma < 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu <= 0)) 
    stop(paste("nu must be positive", "\n", ""))
  
  loglik <- log(nu) + log(mu) + sigma*x - mu/sigma*(exp(sigma*x) - 1) + 
    (nu - 1)*log(1 - exp(- mu/sigma*(exp(sigma*x) - 1)))
  
  if (log == FALSE) 
    density <- exp(loglik) 
  else 
    density <- loglik
  return(density)
}
#' @export
#' @rdname dGGD
pGGD <- function(q, mu, sigma, nu, lower.tail=TRUE, log.p=FALSE){
  if (any(q < 0)) 
    stop(paste("q must be positive", "\n", ""))
  if (any(mu <= 0)) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma < 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu <= 0)) 
    stop(paste("nu must be positive", "\n", ""))
  
  cdf <- (1 - exp((-mu/sigma)*(exp(sigma*q)-1)))^nu
  
  if (lower.tail == TRUE) 
    cdf <- cdf
  else cdf <- 1 - cdf
  if (log.p == FALSE) 
    cdf <- cdf
  else cdf <- log(cdf)
  cdf
}
#' @export
#' @rdname dGGD
qGGD <- function(p, mu, sigma, nu, lower.tail=TRUE, log.p=FALSE){
  if (any(mu <= 0)) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma < 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu <= 0)) 
    stop(paste("nu must be positive", "\n", ""))
  
  if (log.p == TRUE) p <- exp(p)
  else p <- p
  if (lower.tail == TRUE) p <- p
  else p <- 1 - p
  if (any(p < 0) | any(p > 1)) 
    stop(paste("p must be between 0 and 1", "\n", ""))
  
  q <- 1/sigma * log(1 - sigma/mu*log(1 - p^(1/nu)))
  q
}
#' @importFrom stats runif
#' @export
#' @rdname dGGD
rGGD <- function(n, mu, sigma, nu){
  if(any(n <= 0))
    stop(paste("n must be positive","\n",""))
  if (any(mu <= 0)) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma < 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu <= 0)) 
    stop(paste("nu must be positive", "\n", ""))
  
  n <- ceiling(n)
  p <- runif(n)
  r <- qGGD(p, mu, sigma, nu)
  r
}
#' @export
#' @rdname dGGD
hGGD<-function(x, mu, sigma, nu){
  if (any(x < 0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(mu <= 0)) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma < 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu <= 0)) 
    stop(paste("nu must be positive", "\n", ""))
  
  h <- dGGD(x, mu, sigma, nu, log=FALSE) / 
    pGGD(q=x, mu, sigma, nu, lower.tail=FALSE, log.p=FALSE)
  h
}

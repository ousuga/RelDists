#' The Modified Weibull distribution
#' 
#' @author Johan David Marin Benjumea, \email{johand.marin@@udea.edu.co}
#' 
#' @description 
#' Density, distribution function, quantile function, 
#' random generation and hazard function for the modified weibull distribution 
#' with parameters \code{mu}, \code{sigma} and \code{nu}.
#' 
#' @param x,q	vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. 
#' @param mu shape parameter one.    
#' @param sigma parameter two.
#' @param nu scale parameter three.        
#' @param log,log.p	logical; if TRUE, probabilities p are given as log(p).	
#' @param lower.tail logical; if TRUE (default), probabilities are 
#' P[X <= x], otherwise, P[X > x].
#' 
#' @details 
#' The modified weibull distribution with parameters \code{mu}, \code{sigma}
#' and \code{nu} has density given by
#' 
#' \eqn{f(x) = \mu (\sigma + \nu x) x^{\sigma - 1} \exp(\nu x) \exp(-\mu x^{\sigma} \exp(\nu x))}
#' 
#' for \eqn{x > 0}, \eqn{\mu > 0}, \eqn{\sigma \geq 0} and \eqn{\nu \geq 0}.
#' 
#' @return 
#' \code{dMW} gives the density, \code{pMW} gives the distribution 
#' function, \code{qMW} gives the quantile function, \code{rMW}
#' generates random deviates and \code{hMW} gives the hazard function.
#' 
#' @example examples/examples_dMW.R
#' 
#' @references
#'\insertRef{almalki2014modifications}{RelDists}
#'
#'\insertRef{lai2003modified}{RelDists}
#' 
#' @export
dMW <- function(x, mu, sigma, nu, log = FALSE){
  if (any(x < 0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(mu <= 0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma < 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu < 0)) 
    stop(paste("nu must be positive", "\n", ""))
  
  loglik<- log(mu) + log(sigma + nu*x) + (sigma-1)*log(x) +
    nu*x - mu*(x^sigma)*exp(nu*x)
  
  if (log == FALSE) density<- exp(loglik)
  else density <- loglik
  return(density)
}
#' @export
#' @rdname dMW
pMW <- function(q, mu, sigma, nu, lower.tail=TRUE, log.p = FALSE){
  if (any(q < 0)) 
    stop(paste("q must be positive", "\n", ""))
  if (any(mu <= 0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma < 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu < 0)) 
    stop(paste("nu must be positive", "\n", ""))
  
  cdf <- 1- exp(-mu*(q^sigma) * exp(nu*q))
  
  if (lower.tail == TRUE) cdf <- cdf
  else cdf <- 1 - cdf
  if (log.p == FALSE)  cdf <- cdf
  else cdf <- log(cdf)
  cdf
}
#' @export
#' @rdname dMW
qMW <- function(p, mu, sigma, nu, lower.tail = TRUE, log.p = FALSE){
  if (any(mu <= 0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma < 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu < 0)) 
    stop(paste("nu must be positive", "\n", ""))
  
  if (log.p == TRUE) p <- exp(p)
  else p <- p
  if (lower.tail == TRUE) p <- p
  else p <- 1 - p
  if (any(p < 0) | any(p > 1)) 
    stop(paste("p must be between 0 and 1", "\n", ""))
  
  fda <- function(x, mu, sigma, nu){
    1- exp(-mu*(x^sigma) * exp(nu*x))
  }
  fda1 <- function(x, mu, sigma, nu, p) {fda(x, mu,sigma,nu) - p}
  r_de_la_funcion <- function(mu, sigma, nu, p) {
    uniroot(fda1, interval=c(0,99999), mu, sigma, nu, p)$root
  }
  r_de_la_funcion <- Vectorize(r_de_la_funcion)
  q <- r_de_la_funcion(mu, sigma, nu, p)
  q
}
#' @export
#' @rdname dMW
rMW <- function(n, mu, sigma, nu){
  if (any(mu <= 0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma < 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu < 0)) 
    stop(paste("nu must be positive", "\n", ""))
  
  n <- ceiling(n)
  p <- runif(n)
  r <- qMW(p, mu, sigma, nu) 
  r
}
#' @export
#' @rdname dMW
hMW<-function(x, mu, sigma, nu){
  if (any(x < 0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(mu <= 0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma < 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu < 0)) 
    stop(paste("nu must be positive", "\n", ""))
  
  h <- dMW(x, mu, sigma, nu, log = FALSE)/
    pMW(q=x, mu, sigma, nu, lower.tail=FALSE, log.p = FALSE)
  h  
}

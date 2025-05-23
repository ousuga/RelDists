#' The New Exponentiated Exponential distribution
#' 
#' @author Juliana Garcia, \email{juliana.garciav@udea.edu.co}
#' 
#' @description
#' Density, distribution function, quantile function, 
#' random generation and hazard function for the two-parameter 
#' New Exponentiated Exponential with
#' parameters \code{mu} and \code{sigma}.
#' 
#' @param x,q	vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. 
#' @param mu parameter.    
#' @param sigma parameter.
#' @param log,log.p	logical; if TRUE, probabilities p are given as log(p).	
#' @param lower.tail logical; if TRUE (default), probabilities are 
#' \eqn{P[X <= x]}, otherwise, \eqn{P[X > x]}.
#' 
#' @seealso \link{NEE}
#' 
#' @details 
#' The New Exponentiated Exponential distribution with parameters \code{mu} 
#' and \code{sigma} has density given by
#' 
#' \eqn{f(x | \mu, \sigma) = \log(2^\sigma) \mu \exp(-\mu x) (1-\exp(-\mu x))^{\sigma-1} 2^{(1-\exp(-\mu x))^\sigma}, }
#' 
#' for \eqn{x>0}, \eqn{\mu>0} and \eqn{\sigma>0}.
#' 
#' Note: In this implementation we changed the original parameters 
#' \eqn{\theta} for \eqn{\mu} and \eqn{\alpha} for \eqn{\sigma},
#' we did it to implement this distribution within gamlss framework.
#' 
#' @return 
#' \code{dNEE} gives the density, \code{pNEE} gives the distribution 
#' function, \code{qNEE} gives the quantile function, \code{rNEE}
#' generates random deviates and \code{hNEE} gives the hazard function.
#' 
#' @example examples/examples_dNEE.R
#' 
#' @references
#' Hassan, Anwar, I. H. Dar, and M. A. Lone. "A New Class of Probability 
#' Distributions With An Application to Engineering Data." 
#' Pakistan Journal of Statistics and Operation Research 20.2 (2024): 217-231.
#' 
#' @export
dNEE <- function(x, mu=1, sigma=1, log=FALSE) {
  
  if (any(mu <= 0))    stop("parameter mu has to be positive!")
  if (any(sigma <= 0)) stop("parameter sigma has to be positive!")
  
  w <- log(2^sigma)
  r <- 1 - exp(-mu*x)
  res1 <- log(w) + log(mu) - mu*x
  res2 <- (sigma-1) * log(r) + r^sigma * log(2)
  res <- res1 + res2 
  
  if (log == TRUE)
    result <- res
  else
    result <- exp(res)
  return(result)
}
#' @export
#' @rdname dNEE
pNEE <- function(q, mu=1, sigma=1, lower.tail=TRUE, log.p=FALSE){

  if (any(mu <= 0))    stop("parameter mu has to be positive!")
  if (any(sigma <= 0)) stop("parameter sigma has to be positive!")
  
  z <- 1 - exp(-mu*q)
  cdf <- 2^(z^sigma) - 1
  
  if (lower.tail == TRUE) 
    cdf <- cdf
  else 
    cdf = 1 - cdf
  if (log.p == FALSE) 
    cdf <- cdf
  else 
    cdf <- log(cdf)
  cdf <- ifelse(q < 0, 0, cdf)
  return(cdf)
}
#' @export
#' @rdname dNEE
qNEE <- function(p, mu=1, sigma=1, lower.tail=TRUE, log.p=FALSE) {
  
  if (any(mu <= 0))    stop("parameter mu has to be positive!")
  if (any(sigma <= 0)) stop("parameter sigma has to be positive!")
  
  if (log.p == TRUE) {
    p <- exp(p)
  } else {
    p <- p
  }
  if (lower.tail == TRUE) {
    p <- p
  } else {
    p <- 1 - p
  }
  p <- (-1/mu) * log(1 - (log(1 + p)/log(2))^(1/sigma))
  return(p)
}
#' @importFrom stats runif
#' @export
#' @rdname dNEE
rNEE <- function(n=1, mu=1, sigma=1) {
  
  if (any(mu <= 0))    stop("parameter mu has to be positive!")
  if (any(sigma <= 0)) stop("parameter sigma has to be positive!")
  
  u <- runif(n, 0, 1)
  x <- (-1/mu) * log(1 - (log(1 + u)/log(2))^(1/sigma))
  return(x)
}
#' @export
#' @rdname dNEE
hNEE <- function(x, mu, sigma, log=FALSE){
  # Juliana debe construir esta funcion
  return("hola")
}


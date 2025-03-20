#' The Extended Odd Frechet-Nadarajah-Haghighi
#' 
#' @author Helber Santiago Padilla, \email{hspadillar@unal.edu.co}
#' 
#' @description 
#' Density, distribution function, quantile function, 
#' random generation and hazard function for the Extended Odd Fr?chet-Nadarajah-Haghighi distribution
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
#' @references
#' Nasiru, S. (2018). Extended Odd Fréchet‐G Family of Distributions 
#' Journal of Probability and Statistics, 2018(1), 2931326.
#' 
#' @details 
#'  Tthe Extended Odd Frechet-Nadarajah-Haghighi \code{mu}, 
#' \code{sigma}, \code{nu} and \code{tau} has density given by
#' 
#' \eqn{f(x)= \frac{\mu\sigma\nu\tau(1+\nu x)^{\sigma-1}e^{(1-(1+\nu x)^\sigma)}[1-(1-e^{(1-(1+\nu x)^\sigma)})^{\mu}]^{\tau-1}}{(1-e^{(1-(1+\nu x)^{\sigma})})^{\mu\tau+1}} e^{-[(1-e^{(1-(1+\nu x)^\sigma)})^{-\mu}-1]^{\tau}},}
#' 
#' for \eqn{x > 0}, \eqn{\mu > 0}, \eqn{\sigma > 0}, \eqn{\nu > 0} and \eqn{\tau > 0}.
#'  
#' @return 
#' \code{dEOFNH} gives the density, \code{pEOFNH} gives the distribution 
#' function, \code{qEOFNH} gives the quantile function, \code{rEOFNH}
#' generates random numbers and \code{hEOFNH} gives the hazard function.
#' 
#' @example examples/examples_dEOFNH.R
#'
#' @export
dEOFNH <- function(x, mu, sigma, nu, tau, log=FALSE){
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
  
  w <- 1 + nu * x
  z <- 1 - exp(1 - w^sigma)
  
  part1 <- log(mu) + log(sigma) + log(nu) + log(tau)
  part2 <- (sigma-1) * log(w) + log(1-z)
  part3 <- (tau-1) * log(1-z^mu) - (mu*tau+1) * log(z)
  part4 <- -(z^(-mu) - 1)^tau
  
  pdf <- part1 + part2 + part3 + part4
  
  if (log)
    pdf <- pdf
  else
    pdf <- exp(pdf)
  
  return(pdf)
}
#' @export
#' @rdname dEOFNH
pEOFNH <- function(q, mu, sigma, nu, tau, lower.tail=TRUE, log.p=FALSE){
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
  
  term1 <- 1+nu*q
  term2 <- exp(1-(term1)^sigma)
  term3 <- 1-term2
  term5 <- (1-term3^mu)/(term3^mu)
  cdf <- exp(-(term5)^tau)
  if (lower.tail == TRUE) cdf <- cdf
  else cdf <- 1 - cdf 
  if (log.p == FALSE) cdf <- cdf
  else cdf <- log(cdf)
  cdf
}
#' @export
#' @rdname dEOFNH
qEOFNH <- function(p, mu, sigma, nu, tau, lower.tail=TRUE, log.p=FALSE){
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
  
  term1 <- (-log(p))^(1/tau)
  term2 <- (1+term1)^(-1/mu)
  term3 <- log(1-term2)
  term4 <- (1-term3)^(1/sigma)
  term5 <- (term4-1)/(nu)
  q <- term5
  q
}
#' @importFrom stats runif
#' @export
#' @rdname dEOFNH
rEOFNH <- function(n, mu, sigma, nu, tau){
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
  r <- qEOFNH(p, mu, sigma, nu, tau)
  r
}
#' @export
#' @rdname dEOFNH
hEOFNH <- function(x, mu, sigma, nu, tau){
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
  
  h <- dEOFNH(x, mu, sigma, nu, tau, log=FALSE) / pEOFNH(x, mu, sigma, nu, tau, lower.tail=FALSE, log.p=FALSE)
  h
}  

#' The Extended Odd Frechet-Nadarajah-Haghighi
#' 
#' @author Helber Santiago Padilla
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
#' @param lower.tail logical; if TRUE (default), probabilities are P[X <= x], otherwise, P[X > x].
#' 
#' @details 
#'  Tthe Extended Odd Fr?chet-Nadarajah-Haghighi \code{mu}, 
#' \code{sigma}, \code{nu} and \code{tau} has density given by
#' 
#' \eqn{f(x)= \frac{\mu\sigma\nu\tau(1+\nu x)^{\sigma-1}e^{(1-(1+\nu x)^\sigma)}[1-(1-e^{(1-(1+\nu x)^\sigma)})^{\mu}]^{\tau-1}}{(1-e^{(1-(1+\nu x)^{\sigma})})^{\mu\tau+1}} e^{-[(1-e^{(1-(1+\nu x)^\sigma)})^{-\mu}-1]^{\tau}},}
#' 
#' for \eqn{x > 0}, \eqn{\mu > 0}, \eqn{\sigma > 0}, \eqn{\nu > 0} and \eqn{\tau > 0}. 
#' @return 
#' \code{dEOFNH} gives the density, \code{pEOFNH} gives the distribution 
#' function, \code{qEOFNH} gives the quantile function, \code{rEOFNH}
#' generates random numbers and \code{hEOFNH} gives the hazard function.
#' 
#' @examples 
#' ##The probability density function
#' par(mfrow=c(1,1))
#'  curve(dEOFNH(x, mu=18.5, sigma=5.1, nu=0.1, tau=0.1), from=0, to=10,
#'      ylim=c(0, 0.25), col="red", las=1, ylab="f(x)")
#' 
#' ## The cumulative distribution and the Reliability function
#' par(mfrow = c(1, 2))
#' curve(pEOFNH(x,mu=18.5, sigma=5.1, nu=0.1, tau=0.1), from = 0, to = 10, 
#' ylim = c(0, 1), col = "red", las = 1, ylab = "The cumulative distribution function")
#' curve(pEOFNH(x, mu=18.5, sigma=5.1, nu=0.1, tau=0.1, lower.tail = FALSE), 
#' from = 0, to = 10, ylim = c(0, 1), col = "red", las = 1, ylab = "The Reliability function")
#' 
#' ##The quantile function
#' p <- seq(from=0, to=0.99999, length.out=100)
#' plot(x=qEOFNH(p, mu=18.5, sigma=5.1, nu=0.1, tau=0.1), y=p, xlab="Quantile",
#'      las=1, ylab="Probability")
#' curve(pEOFNH(x, mu=18.5, sigma=5.1, nu=0.1, tau=0.1), from=0, add=TRUE, col="red")
#' 
#' ##The random function
#' hist(rEOFNH(n=10000, mu=18.5, sigma=5.1, nu=0.1, tau=0.1), freq=FALSE,
#'      xlab="x", las=1, main="")
#' curve(dEOFNH(x, mu=18.5, sigma=5.1, nu=0.1, tau=0.1), from=0, add=TRUE, col="red", ylim=c(0,1.25))
#' 
#' ##The Hazard function
#' par(mfrow=c(1,1))
#' curve(hEOFNH(x, mu=18.5, sigma=5.1, nu=0.1, tau=0.1), from=0, to=10, ylim=c(0, 1),
#'      col="red", ylab="Hazard function", las=1)
#'      
#' @references
#'\insertRef{nasiru2018extended}{RelDists}
#'
#' @export

dEOFNH <- function(x, mu, sigma, nu, tau){
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
  term1 <- 1+nu*x
  term2 <- exp(1-(term1)^sigma)
  term3 <- 1-term2
  term4 <- exp(-(term3^(-mu)-1)^tau)
  density <- term4*(mu*sigma*tau*nu*term1^(sigma-1)*term2*(1-term3^mu)^(tau-1))/(term3^(mu*tau+1))
  return(density)
}
#' @export
#' @rdname dEOFNH

pEOFNH <- function(x, mu, sigma, nu, tau, lower.tail=TRUE){
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
  
  term1 <- 1+nu*x
  term2 <- exp(1-(term1)^sigma)
  term3 <- 1-term2
  term5 <- (1-term3^mu)/(term3^mu)
  cdf <- exp(-(term5)^tau)
  if (lower.tail == TRUE) cdf <- cdf
  else cdf <- 1 - cdf 
  cdf
}

#' @export
#' @rdname dEOFNH


qEOFNH <- function(p, mu, sigma, nu, tau, lower.tail=TRUE){
  if (any(mu <= 0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu <= 0)) 
    stop(paste("nu must be positive", "\n", ""))
  if (any(tau <= 0)) 
    stop(paste("tau must be positive", "\n", ""))
  if (any(p < 0) | any(p > 1)) 
    stop(paste("p must be between 0 and 1", "\n", ""))
  if (lower.tail == TRUE) p <- p
  else p <- 1 - p 
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
  
  h <- dEOFNH(x, mu, sigma, nu, tau) / pEOFNH(x, mu, sigma, nu, tau, lower.tail=FALSE)
  h
}  

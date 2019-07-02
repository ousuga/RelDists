#' The Sarhan and Zaindin's Modified Weibull distribution
#' 
#' @author Johan David Marin Benjumea, \email{johand.marin@@udea.edu.co}
#' 
#' @description 
#' Density, distribution function, quantile function, 
#' random generation and hazard function for Sarhan and Zaindins modified weibull distribution
#' with parameters \code{mu}, \code{sigma} and \code{nu}.
#' 
#' @param x,q	vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. 
#' @param mu scale parameter.    
#' @param sigma shape parameter.
#' @param nu shape parameter.
#' @param log,log.p	logical; if TRUE, probabilities p are given as log(p).	
#' @param lower.tail logical; if TRUE (default), probabilities are 
#' P[X <= x], otherwise, P[X > x].
#' @details 
#' The Sarhan and Zaindins modified weibull with parameters \code{mu}, 
#' \code{sigma} and \code{nu} has density given by
#' 
#' \eqn{f(x)=(\mu + \sigma \nu x^(\nu - 1)) \exp(- \mu x - \sigma x^\nu)}
#' 
#' for \eqn{x > 0}, \eqn{\mu > 0}, \eqn{\sigma > 0} and \eqn{\nu > 0}. 
#' 
#' @return 
#' \code{dSZMW} gives the density, \code{pSZMW} gives the distribution 
#' function, \code{qSZMW} gives the quantile function, \code{rSZMW}
#' generates random deviates and \code{hSZMW} gives the hazard function.
#' 
#' @examples 
#' 
#' ## The probability density function
#' curve(dSZMW(x, mu = 2, sigma = 1.5, nu = 0.2), from = 0, to = 2, 
#'       ylim = c(0, 1.7), col = "red", las = 1, ylab = "The probability density function")
#' ## The cumulative distribution and the Reliability function
#' par(mfrow = c(1, 2))
#' curve(pSZMW(x, mu = 2, sigma = 1.5, nu = 0.2), from = 0, to = 2, ylim = c(0, 1),
#'       col = "red", las = 1, ylab = "The cumulative distribution function")
#' curve(pSZMW(x, mu = 2, sigma = 1.5, nu = 0.2, lower.tail = FALSE), from = 0,
#'       to = 2, ylim = c(0, 1), col = "red", las = 1, ylab = "The Reliability function")
#' 
#' ## The quantile function
#' p <- seq(from = 0, to = 0.99999, length.out = 100)
#' plot(x = qSZMW(p = p, mu = 2, sigma = 1.5, nu = 0.2), y = p, xlab = "Quantile",
#'      las = 1, ylab = "Probability")
#' curve(pSZMW(x, mu = 2, sigma = 1.5, nu = 0.2), from = 0, add = TRUE, col = "red")
#' 
#' ## The random function
#' hist(rSZMW(n = 1000, mu = 2, sigma = 1.5, nu = 0.2), freq = FALSE, xlab = "x",
#'      las = 1, main = "")
#' curve(dSZMW(x, mu = 2, sigma = 1.5, nu = 0.2),  from = 0, add = TRUE, col = "red")
#' 
#' ## The Hazard function
#' par(mfrow=c(1,1))
#' curve(hSZMW(x, mu = 2, sigma = 1.5, nu = 0.2), from = 0, to = 3, ylim = c(0, 8),
#'       col = "red", ylab = "The hazard function", las = 1)
#' 
#' @export
dSZMW<-function(x, mu, sigma, nu, log=FALSE){
  if (any(x < 0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(mu <= 0)) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu <= 0)) 
    stop(paste("nu must be positive", "\n", ""))
  
  loglik<- log(mu + sigma*nu*x^(nu-1)) - mu*x - sigma*x^nu
  
  if (log == FALSE) density<- exp(loglik)
  else density <- loglik
  return(density)
}
#' @export
#' @rdname dSZMW
pSZMW <- function(q, mu, sigma, nu, lower.tail=TRUE, log.p=FALSE){
  if (any(q < 0)) 
    stop(paste("q must be positive", "\n", ""))
  if (any(mu <= 0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu <= 0)) 
    stop(paste("nu must be positive", "\n", ""))
  
  cdf <- 1- exp(-mu*q -sigma*(q^nu))
  if (lower.tail == TRUE) cdf <- cdf
  else cdf <- 1 - cdf
  if (log.p == FALSE) cdf <- cdf
  else cdf <- log(cdf)
  cdf
}
#' @export
#' @rdname dSZMW
qSZMW <- function(p, mu, sigma, nu, lower.tail=TRUE, log.p=FALSE){
  if (any(mu <= 0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu <=0 )) 
    stop(paste("nu must be positive", "\n", ""))
  
  if (log.p == TRUE) p <- exp(p)
  else p <- p
  if (lower.tail == TRUE) p <- p
  else p <- 1 - p
  if (any(p < 0) | any(p > 1)) 
    stop(paste("p must be between 0 and 1", "\n", ""))
  
  fda <- function(x, mu, sigma, nu){ 1- exp(-mu*x - sigma*(x^nu))}
  fda1 <- function(x, mu, sigma, nu, p) {fda(x, mu, sigma, nu) - p}
  r_de_la_funcion <- function(mu, sigma, nu, p) {
    uniroot(fda1, interval=c(0,1e+06), mu, sigma, nu, p)$root
  }
  r_de_la_funcion <- Vectorize(r_de_la_funcion)
  q <- r_de_la_funcion(mu,sigma,nu, p)
  q
}
#' @export
#' @rdname dSZMW
rSZMW<- function(n, mu, sigma, nu){
  if (any(mu <= 0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu <= 0)) 
    stop(paste("nu must be positive", "\n", ""))
  
  n <- ceiling(n)
  p <- runif(n)
  r <- qSZMW(p, mu, sigma, nu)
  r
}
#' @export
#' @rdname dSZMW
hSZMW<-function(x, mu, sigma, nu){
  if (any(x < 0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(mu <= 0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu <= 0)) 
    stop(paste("nu must be positive", "\n", ""))
  
  h <- dSZMW(x, mu, sigma, nu, log=FALSE)/
    pSZMW(q=x, mu, sigma, nu, lower.tail=FALSE, log.p=FALSE)
  h  
}
#' The Exponentiated Exponential Binomial distribution
#' 
#' @author Jaime Mosquera Gutiérrez, \email{jmosquerag@@unal.edu.co}
#' 
#' @description 
#' Density, distribution function, quantile function, 
#' random generation and hazard function for the exponentiated exponential binomial distribution with
#' parameters \code{mu}, \code{sigma} and \code{nu}.
#' 
#' @param x,q	vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. 
#' @param mu scale parameter.
#' @param sigma,nu shape parameters.
#' @param size number of trials (one or more).
#' @param log,log.p	logical; if TRUE, probabilities p are given as log(p).	
#' @param lower.tail logical; if TRUE (default), probabilities are P[X <= x], otherwise, P[X > x].
#' 
#' @details 
#' The exponentiated exponential binomial distribution with parameters \code{size =}\eqn{m}, 
#' \code{mu}, \code{sigma} and \code{nu} has density given by
#' 
#' \eqn{f(x) = \sigma\nu\mu m \exp(-\sigma x) (1 - \exp(-\sigma x) )^(\nu-1) (1 - \mu(1 - \exp(-\sigma x))^\nu)^(m-1) ,}
#' 
#' for \eqn{x > 0}, \eqn{0 < \mu < 1}, \eqn{\sigma > 0} and \eqn{\nu >0}.
#' 
#' @return 
#' \code{dEEB} gives the density, \code{pEEB} gives the distribution 
#' function, \code{qEEB} gives the quantile function, \code{rEEB}
#' generates random deviates and \code{hEEB} gives the hazard function.
#'
#' @examples  
#' ## The probability density function
#' par(mfrow = c(1, 2))
#' curve(dEEB(x, mu=1, sigma=1, nu=0.5, size=10), from=0, to=1.5, ylim=c(0, 1.5), 
#'       col="red", las=1, ylab="The probability density function")
#' curve(dEEB(x, mu=1, sigma=1, nu=5, size=10), from=0, to=1.5, ylim=c(0, 3), 
#'       col="red", las=1, ylab="The probability density function")
#'       
#' ## The cumulative distribution and the Reliability function
#' par(mfrow = c(1, 2))
#' curve(pEEB(x, mu=1, sigma=1, nu=5, size=10), from=0, to=4, ylim=c(0, 1), 
#'       col="red", las = 1, ylab = "The cumulative distribution function")
#' curve(pEEB(x, mu=1, sigma=1, nu=5, size=10, lower.tail=FALSE), from=0, to=4,  
#'       ylim=c(0, 1), col="red", las=1, ylab="The Reliability function")
#' 
#' ## The quantile function
#' p <- seq(from=0, to=0.998, length.out=100)
#' plot(x=qEEB(p, mu=1, sigma=1, nu=5, size=10), y=p, xlab="Quantile", 
#'     las=1, ylab="Probability")
#' curve(pEEB(x, mu=1, sigma=1, nu=5, size=10), from=0, add=TRUE, col="red")
#' 
#' ## The random function
#' hist(rEEB(n=10000, mu=1, sigma=1, nu=5, size=10), freq=FALSE, 
#'      ylim = c(0,3),xlab="x", las=1, main="")
#' curve(dEEB(x, mu=1, sigma=1, nu=5, size=10), from=0, ylim=c(0, 3), 
#'       add=TRUE, col="red")
#' 
#' ## The Hazard function
#' par(mfrow=c(1,2))
#' curve(hEEB(x, mu=1, sigma=1, nu=0.5, size=10), from=0, to=4, ylim=c(0, 3), 
#'       col="red", ylab="The hazard function", las=1)
#' curve(hEEB(x, mu=1, sigma=1, nu=5, size=10), from=0, to=6, ylim=c(0, 3), 
#'       col="red", ylab="The hazard function", las=1)
#'       
#' @references
#' \insertRef{Bakouch2012}{RelDists}
#'
#' @importFrom Rdpack reprompt
#'
#' @export
dEEB <- function(x, mu, sigma, nu, size, log=FALSE){
  if (any(x<0))
    stop(paste("x must greater than zero", "\n", ""))
  if (any(mu < 0) | any(mu > 1))
    stop(paste("mu must be between 0 and 1", "\n", ""))
  if (any(sigma<=0 ))
    stop(paste("sigma must be greater than zero", "\n", ""))
  if (any(nu<=0 ))
    stop(paste("nu must be greater than zero", "\n", ""))
  
  loglik <- log(sigma * nu * size *  mu) - sigma * x +
    (nu - 1) * log(1 - exp(-sigma * x)) -
    log(1 - (1 - mu)^size) +
    (size - 1) * log(1 - mu * (1 - exp(-sigma * x))^nu)
  if (log == FALSE)
    density <- exp(loglik) else density <- loglik
  return(density)
}
#' @export
#' @rdname dEEB
pEEB <- function(q, mu, sigma, nu, size, lower.tail=TRUE, log.p=FALSE){
  if (any(q < 0)) 
    stop(paste("q must be positive", "\n", ""))
  if (any(mu < 0) | any(mu > 1))
    stop(paste("mu must be between 0 and 1", "\n", ""))
  if (any(sigma<=0 ))
    stop(paste("sigma must be greater than zero", "\n", ""))
  if (any(nu<=0 ))
    stop(paste("nu must be greater than zero", "\n", ""))

  cdf <- ( 1 - (1 - mu*(1 - exp(-sigma*q))^nu)^size )/( 1 - (1 - mu)^size )
  
  if (lower.tail == TRUE) cdf <- cdf
  else cdf <- 1 - cdf
  if (log.p == FALSE) cdf <- cdf
  else cdf <- log(cdf)
  cdf
}
#' @export
#' @rdname dEEB
qEEB <-  function(p, mu, sigma, nu, size, lower.tail=TRUE, log.p=FALSE){
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
  else p <- 1 - p
  if (any(p < 0) | any(p > 1))
    stop(paste("p must be between 0 and 1", "\n", ""))
  
  q <- -(1/sigma) * log( 1 - ( (1/mu)*(1 - (1 - (1 - (1 - mu)^size)*p)^(1/size)) )^(1/nu) )
  q
}
#' @export
#' @rdname dEEB
rEEB <- function(n, mu, sigma, nu, size){
  if(any(n <= 0))
    stop(paste("n must be a positive integer","\n",""))
  if (any(mu <= 0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu <= 0)) 
    stop(paste("nu must be positive", "\n", ""))
  
  n <- ceiling(n)
  p <- runif(n)
  r <- qEEB(p, mu, sigma, nu, size)
  r
}
#' @export
#' @rdname dEEB
hEEB <- function(x, mu, sigma, nu, size){
  if (any(x < 0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(mu <= 0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma*nu <= 0)) 
    stop(paste("Product sigma*nu must be positive", "\n", ""))
  
  h <- dEEB(x, mu, sigma, nu, size, log=FALSE)/
    pEEB(q=x, mu, sigma, nu, size, lower.tail=FALSE, log.p=FALSE)
  h
}
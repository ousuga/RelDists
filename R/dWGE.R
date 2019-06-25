#' The Weighted Gamma-Exponential distribution
#' 
#' @author Jaime Mosquera Gutiérrez, \email{jmosquerag@@unal.edu.co}
#' 
#' @description 
#' Density, distribution function, quantile function, 
#' random generation and hazard function for the Weighted Gamma-Exponential 
#' distribution with parameters \code{mu}, \code{sigma} and \code{nu}.
#' 
#' @param x,q	vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. 
#' @param mu scale parameter.
#' @param sigma,nu shape parameters.
#' @param log,log.p	logical; if TRUE, probabilities p are given as log(p).	
#' @param lower.tail logical; if TRUE (default), probabilities are P[X <= x], otherwise, P[X > x].
#' 
#' @details 
#' The Weighted Gamma-Exponential distribution with parameters \code{mu},
#' \code{sigma} and \code{nu} has density given by
#' 
#' \eqn{f(x)= \mu^\nu x^(\nu-1) \exp(-\mu x)(1 - \exp(-\sigma\mu x)) / \Gamma(\nu)(1 - (1 + \alpha)^(-\nu) ,}
#' 
#' for \eqn{x > 0}, \eqn{\mu > 0}, \eqn{\sigma > 0} and \eqn{\nu >0}.
#' 
#' @return 
#' \code{dWGE} gives the density, \code{pWGE} gives the distribution 
#' function, \code{qWGE} gives the quantile function, \code{rdWGE}
#' generates random deviates and \code{hWGE} gives the hazard function.
#' 
#' @examples  
#' ## The probability density function 
#' par(mfrow = c(1,2))
#' curve(dWGE(x, mu=1, sigma=1, nu=1), from=0, to=8, ylim=c(0, 0.6), 
#'       col="red", las=1, ylab="The probability density function")
#' curve(dWGE(x, mu=1, sigma=1, nu=0.5), from=0, to=8, ylim=c(0, 0.8), 
#'       col="red", las=1, ylab="The probability density function")
#' 
#' ## The cumulative distribution and the Reliability function
#' par(mfrow = c(1, 2))
#' curve(pWGE(x, mu=1, sigma=1, nu=1), from=0, to=10, ylim=c(0, 1), 
#'       col="red", las = 1, ylab = "The cumulative distribution function")
#' curve(pWGE(x, mu=1, sigma=1, nu=1, lower.tail=FALSE), from=0, to=10,  
#'       ylim=c(0, 1), col="red", las=1, ylab="The Reliability function")
#'
#' ## The quantile function
#' p <- seq(from=0, to=0.998, length.out=100)
#' plot(x=qWGE(p, mu=1, sigma=1, nu=1), y=p, xlab="Quantile", 
#'     las=1, ylab="Probability")
#' curve(pWGE(x, mu=1, sigma=1, nu=1), from=0, add=TRUE, col="red")
#' 
#' ## The random function
#' hist(rWGE(n=100, mu=1, sigma=1, nu=1), freq=FALSE, 
#'      ylim = c(0,0.6),xlab="x", las=1, main="")
#' curve(dWGE(x, mu=1, sigma=1, nu=1), from=0, ylim=c(0, 0.6), 
#'       add=TRUE, col="red")
#' 
#' ## The Hazard function
#' par(mfrow=c(1,2))
#' curve(hWGE(x, mu=1, sigma=1, nu=1), from=0, to=8, ylim=c(0, 1.3), 
#'       col="red", ylab="The hazard function", las=1)
#' curve(hWGE(x, mu=1, sigma=10, nu=0.5), from=0.1, to=8, ylim=c(0, 1.3), 
#'       col="red", ylab="The hazard function", las=1)
#'
#' @export
dWGE <- function(x, mu, sigma, nu, log=FALSE){
  if (any(x < 0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(mu <= 0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu <= 0)) 
    stop(paste("nu must be positive", "\n", ""))
  
  loglik <- nu * log(nu) + (nu - 1) * log(x) - nu * x - log(base::gamma(nu)) +
    log(1 - exp(-sigma*mu*x)) - log(1 - (1 + sigma)^(-nu))
  
  if (log == FALSE)
    density<- exp(loglik)
  else 
    density <- loglik
  return(density)
}
#' @export
#' @rdname dWGE
pWGE <- function(q, mu, sigma, nu, lower.tail=TRUE, log.p=FALSE){
  if (any(q < 0))
    stop(paste("q must be positive", "\n", ""))
  if (any(mu <= 0 ))
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0))
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu <= 0))
    stop(paste("nu must be positive", "\n", ""))

  # The incomplete gamma function
  igamma <- function(a, x) {
    stats::pgamma(x, shape=a, scale=1, lower.tail=TRUE) * base::gamma(a)
  }

  cdf <- ( (1 + sigma)^nu * igamma(a=nu, x=mu*q) -
           igamma(a=nu, x=mu*(1 + sigma)*q) ) / ( gamma(nu)*((1 + sigma)^nu - 1) )

  if (lower.tail == TRUE)
    cdf <- cdf
  else cdf <- 1 - cdf
  if (log.p == FALSE)
    cdf <- cdf
  else cdf <- log(cdf)
  cdf
}
#' @export
#' @rdname dWGE
qWGE <- function(p, mu, sigma, nu, lower.tail=TRUE, log.p=FALSE){
  if (any(mu <= 0 ))
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0))
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu <= 0))
    stop(paste("nu must be positive", "\n", ""))

  if (log.p)
    p <- exp(p)
  else p <- p
  if (lower.tail)
    p <- p
  else p <- 1 - p
  if (any(p < 0) | any(p > 1))
    stop(paste("p must be between 0 and 1", "\n", ""))

  F.inv <- function(y, mu, sigma, nu) {
    uniroot(function(x) {pWGE(x, mu, sigma, nu) - y},
            interval=c(0, 99999))$root
  }
  F.inv <- Vectorize(F.inv)
  F.inv(p, mu, sigma, nu)
}
#' @importFrom stats runif
#' @export
#' @rdname dWGE
rWGE <- function(n, mu, sigma, nu){
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
  r <- qWGE(p, mu, sigma, nu)
  r
}
#' @export
#' @rdname dWGE
hWGE <- function(x, mu, sigma, nu){
  if (any(x < 0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(mu <= 0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma*nu <= 0)) 
    stop(paste("Product sigma*nu must be positive", "\n", ""))
  
  h <- dWGE(x, mu, sigma, nu, log=FALSE)/
       pWGE(q=x, mu, sigma, nu, lower.tail=FALSE, log.p=FALSE)
  h
}


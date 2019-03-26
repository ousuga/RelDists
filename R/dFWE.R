#' The Flexible Weibull Extension distribution
#' 
#' @description
#' Density, distribution function, quantile function, 
#' random generation and hazard function for the Flexible Weibull Extension distribution with
#' parameters \code{mu} and \code{sigma}.
#' 
#' @param x,q	vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. 
#' @param mu parameter.    
#' @param sigma parameter.
#' @param log,log.p	logical; if TRUE, probabilities p are given as log(p).	
#' @param lower.tail logical; if TRUE (default), probabilities are P[X <= x], otherwise, P[X > x].
#' 
#' @details 
#' The Flexible Weibull extension with parameters \code{mu} and \code{sigma}
#' has density given by
#' 
#' \eqn{f(x) = (\mu + \sigma/x^2) \exp(\mu x - \sigma/x) \exp(-\exp(\mu x-\sigma/x))}
#' 
#' for x>0.
#' 
#' @return 
#' \code{dFWE} gives the density, \code{pFWE} gives the distribution 
#' function, \code{qFWE} gives the quantile function, \code{rFWE}
#' generates random deviates and \code{hFWE} gives the hazard function.
#' 
#' @examples  
#' ## The probability density function
#' curve(dFWE(x, mu=0.75, sigma=0.5), from=0, to=3, 
#'       ylim=c(0, 1.7), col="red", las=1, ylab="f(x)")
#' 
#' ## The cumulative distribution and the Reliability function
#' par(mfrow=c(1, 2))
#' curve(pFWE(x, mu=0.75, sigma=0.5), from=0, to=3, 
#'       col="red", las=1, ylab="F(x)")
#' curve(pFWE(x, mu=0.75, sigma=0.5, lower.tail=FALSE), 
#'       from=0, to=3, col="red", las=1, ylab="S(x)")
#' 
#' ## The quantile function
#' p <- seq(from=0, to=0.99999, length.out=100)
#' plot(x=qFWE(p, mu=0.75, sigma=0.5), y=p, xlab="Quantile",
#'      las=1, ylab="Probability")
#' curve(pFWE(x, mu=0.75, sigma=0.5), from=0, add=TRUE, col="red")
#' 
#' ## The random function
#' hist(rFWE(n=1000, mu=2, sigma=0.5), freq=FALSE, xlab="x", 
#'      ylim=c(0, 2), las=1, main="")
#' curve(dFWE(x, mu=2, sigma=0.5), from=0, to=3, add=TRUE, col="red")
#' 
#' ## The Hazard function
#' par(mfrow=c(1,1))
#' curve(hFWE(x, mu=0.75, sigma=0.5), from=0, to=2, ylim=c(0, 2.5), 
#'       col="red", ylab="Hazard function", las=1)
#' 
#' @export
dFWE <- function(x, mu, sigma, log=FALSE){
  if (any(x < 0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(mu <= 0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  
  loglik <- log(mu + (sigma/x^2)) + (mu*x) - (sigma/x) - 
    exp(mu*x - (sigma/x))
  
  if (log == FALSE) 
    density <- exp(loglik)
  else density <- loglik
  return(density)
}
#' @export
#' @rdname dFWE
pFWE <- function(q, mu, sigma, lower.tail=TRUE, log.p=FALSE){
  if (any(q < 0)) 
    stop(paste("q must be positive", "\n", ""))
  if (any(mu <= 0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  
  cdf <- 1- exp(-exp(mu*q - sigma/q))
  
  if (lower.tail == TRUE) 
    cdf <- cdf
  else cdf <- 1 - cdf
  if (log.p == FALSE) 
    cdf <- cdf
  else cdf <- log(cdf)
  cdf
}
#' @importFrom stats uniroot
#' @export
#' @rdname dFWE
qFWE <- function(p, mu, sigma, lower.tail=TRUE, log.p=FALSE) {
  if (any(mu <= 0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  
  if (log.p == TRUE) 
    p <- exp(p)
  else p <- p
  if (lower.tail == TRUE) 
    p <- p
  else p <- 1 - p
  if (any(p < 0) | any(p > 1)) 
    stop(paste("p must be between 0 and 1", "\n", ""))
  
  fda <- function(x,mu, sigma){
    1- exp(-exp(mu*x - sigma/x))
  }
  fda1 <- function(x, mu, sigma, p) {fda(x, mu, sigma) - p}
  r_de_la_funcion <- function(mu, sigma, p) {
    uniroot(fda1, interval=c(0, 1e+06), mu, sigma, p)$root
  }
  r_de_la_funcion <- Vectorize(r_de_la_funcion)
  q <- r_de_la_funcion(mu, sigma, p)
  q
}
#' @export
#' @rdname dFWE
rFWE <- function(n, mu, sigma){
  if (any(mu <= 0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  
  n <- ceiling(n)
  p <- runif(n)
  r <- qFWE(p, mu, sigma)
  r
}
#' @export
#' @rdname dFWE
hFWE <- function(x, mu, sigma){
  if (any(x < 0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(mu <= 0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  
  h <- dFWE(x, mu, sigma, log = FALSE) / 
    pFWE(q=x, mu, sigma, lower.tail=FALSE, log.p=FALSE)
  h
}
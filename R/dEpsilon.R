#' The Epsilon distribution
#' 
#' @author Daniel Betancur
#' 
#' @description 
#' Density, distribution function, quantile function, 
#' random generation and hazard function for the Epsilon distribution 
#' with parameters \code{mu} and \code{sigma}.
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
#' The Epsilon Distribution with parameters \code{mu} and 
#' \code{sigma} has density given by
#' 
#' \eqn{f(x)= \mu \frac{\sigma^2}{\sigma^2-x^2}\varepsilon_{-\mu,\sigma}\left(x\right),}
#' 
#' where
#' 
#' \eqn{\varepsilon_{-\mu,\sigma}(x) = \left(\frac{x+\sigma}{\sigma-x}\right)^{-\mu\frac{\sigma}{2}},} 
#' 
#' for \eqn{0 < x < \sigma}, \eqn{\mu > 0}  and \eqn{\sigma > 0}. 
#' 
#' @return 
#' \code{dEpsilon} gives the density, \code{pEpsilon} gives the distribution 
#' function, \code{qEpsilon} gives the quantile function, \code{rEpsilon}
#' generates random deviates and \code{hEpsilon} gives the hazard function.
#'
#' @examples  
#' ## The probability density function
#' par(mfrow=c(1,1))
#' curve(dEpsilon(x, mu=0.7, sigma=5), from=0, to=4.999,
#'       ylim=c(0, 0.7), col="red", las=1, ylab="f(x)")
#' 
#' ## The cumulative distribution and the Reliability function
#' par(mfrow=c(1, 2))
#' curve(pEpsilon(x, mu=0.7, sigma=5),
#'       from=0, to=4.999,  col="red", las=1, ylab="F(x)")
#' curve(pEpsilon(x, mu=0.7, sigma=5, lower.tail=FALSE),
#'       from=0, to=4.999, col="red", las=1, ylab="S(x)")
#' 
#' ## The quantile function
#' p <- seq(from=0, to=0.99999, length.out=100)
#' plot(x=qEpsilon(p, mu=0.7, sigma=5), y=p, xlab="Quantile",
#'      las=1, ylab="Probability")
#' curve(pEpsilon(x, mu=0.7, sigma=10), from=0, add=TRUE, col="red")
#' 
#' ## The random function
#' hist(rEpsilon(n=10000, mu=0.7, sigma=5), freq=FALSE,
#'      xlab="x", las=1, main="")
#' curve(dEpsilon(x, mu=0.7, sigma=10), from=0, add=TRUE, col="red")
#' 
#' ## The Hazard function
#' par(mfrow=c(1,1))
#' curve(hEpsilon(x, mu=0.7, sigma=5), from=0, to=4.999, ylim=c(0, 15),
#'       col="red", ylab="Hazard function", las=1)
#'
#' @references
#'\insertRef{dombi2018epsilon}{RelDists}
#'
#' @export
dEpsilon <- function(x, mu, sigma,log=FALSE){
  if (any(x < 0)) 
    stop(paste("x must be between 0 and sigma", "\n", ""))
  if (any(mu <= 0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  
  Epsilon <- ((x + sigma)/(sigma - x))^(-mu*sigma/2)
  density <- mu*((sigma^2)/(sigma^2-x^2))*Epsilon
  
  if (log == FALSE) 
    density <- density
  else 
    density <- log(density)
  return(density) 
}
#' @export
#' @rdname dEpsilon
pEpsilon <- function(q, mu, sigma, 
                     lower.tail=TRUE, log.p=FALSE){
  if (any(q < 0)) 
    stop(paste("q must be between 0 and sigma", "\n", ""))
  if (any(mu <= 0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  
  Epsilon <- ((q + sigma)/(sigma - q))^(-mu*sigma/2)
  cdf <- 1 - Epsilon
  
  if (lower.tail == TRUE) cdf <- cdf
  else cdf <- 1 - cdf 
  if (log.p == FALSE) cdf <- cdf
  else cdf <- log(cdf)
  cdf
}

#' @export
#' @rdname dEpsilon
qEpsilon <- function(p, mu, sigma, 
                     lower.tail=TRUE, log.p=FALSE){
  if (any(mu <= 0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  
  
  if (log.p == TRUE) p <- exp(p)
  else p <- p
  if (lower.tail == TRUE) p <- 1- p
  else p <- p
  if (any(p < 0) | any(p > 1)) 
    stop(paste("p must be between 0 and 1", "\n", ""))
  
  q <- sigma*(((1-p)^(-2/(mu*sigma))+ 1)/((1-p)^(-2/(mu*sigma)- 1)))
  q
}
#' @importFrom stats runif
#' @export
#' @rdname dEpsilon
rEpsilon <- function(n, mu, sigma){
  if(any(n <= 0))
    stop(paste("n must be positive","\n",""))
  if (any(mu <= 0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  
  n <- ceiling(n)
  p <- runif(n)
  r <- qEpsilon(p, mu, sigma)
  r
}
#' @export
#' @rdname dEpsilon
hEpsilon <- function(x, mu, sigma){
  if (any(x < 0)) 
    stop(paste("x must be between 0 and sigma", "\n", ""))
  if (any(mu <= 0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  
  
  h <- dEpsilon(x, mu, sigma,log=FALSE) / 
    pEpsilon(x, mu, sigma,lower.tail=FALSE, log.p=FALSE)
  h
}
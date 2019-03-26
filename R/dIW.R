#' The Inverse Weibull distribution
#' 
#' @description 
#' Density, distribution function, quantile function, 
#' random generation and hazard function for the inverse weibull distribution with
#' parameters \code{mu} and \code{sigma}.
#' 
#' @param x,q	vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. 
#' @param mu scale parameter.
#' @param sigma shape parameters.
#' @param log,log.p	logical; if TRUE, probabilities p are given as log(p).	
#' @param lower.tail logical; if TRUE (default), probabilities are P[X <= x], otherwise, P[X > x].
#' 
#' @details 
#' The inverse weibull distribution with parameters \code{mu} and
#' \code{sigma} has density given by
#' 
#' \eqn{f(x) = \mu \sigma x^{-\sigma-1} \exp(\mu x^{-\sigma})}
#' 
#' for x > 0. 
#' 
#' @return 
#' \code{dIW} gives the density, \code{pIW} gives the distribution 
#' function, \code{qIW} gives the quantile function, \code{rIW}
#' generates random deviates and \code{hIW} gives the hazard function.
#'
#' @examples  
#' ## The probability density function
#' curve(dIW(x, mu=5, sigma=2.5), from=0, to=10,
#'       ylim=c(0, 0.55), col="red", las=1, ylab="f(x)")
#'
#' ## The cumulative distribution and the Reliability function
#' par(mfrow=c(1, 2))
#' curve(pIW(x, mu=5, sigma=2.5),
#'       from=0, to=10,  col="red", las=1, ylab="F(x)")
#' curve(pIW(x, mu=5, sigma=2.5, lower.tail=FALSE),
#'       from=0, to=10, col="red", las=1, ylab="S(x)")
#'             
#' ## The quantile function
#' p <- seq(from=0, to=0.99999, length.out=100)
#' plot(x=qIW(p, mu=5, sigma=2.5), y=p, xlab="Quantile",
#'   las=1, ylab="Probability")
#' curve(pIW(x, mu=5, sigma=2.5), from=0, add=TRUE, col="red")
#'   
#' ## The random function
#' hist(rIW(n=10000, mu=5, sigma=2.5), freq=FALSE, xlim=c(0,60),
#'   xlab="x", las=1, main="")
#' curve(dIW(x, mu=5, sigma=2.5), from=0, add=TRUE, col="red")
#' 
#' ## The Hazard function
#' par(mfrow=c(1,1))
#' curve(hIW(x, mu=5, sigma=2.5), from=0, to=15, ylim=c(0, 0.9),
#'    col="red", ylab="Hazard function", las=1)
#'
#' @export
dIW <- function(x, mu, sigma, log=FALSE){
  if (any(x < 0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(mu <= 0)) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  
  loglik <- log(mu*sigma) - (sigma+1)*log(x) - mu*(x^-sigma)
  
  if (log == FALSE) 
    density <- exp(loglik) 
  else density <- loglik
  return(density)  
}
#' @export
#' @rdname dIW
pIW <- function(q, mu, sigma, lower.tail=TRUE, log.p=FALSE){
  if (any(q < 0)) 
    stop(paste("q must be positive", "\n", ""))
  if (any(mu <= 0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  
  cdf <- exp((-mu)*(q^(-sigma)))
  if (lower.tail == TRUE) 
    cdf <- cdf
  else cdf <- 1 - cdf
  if (log.p == FALSE) 
    cdf <- cdf
  else cdf <- log(cdf)
  cdf
}
#' @export
#' @rdname dIW
qIW <- function(p, mu, sigma, lower.tail = TRUE, log.p = FALSE){
  if (any(mu <= 0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  
  if (log.p == TRUE) 
    p <- exp(p)
  else p <- p
  if (lower.tail == TRUE) 
    p <- p
  else  p <- 1 - p
  if (any(p < 0) | any(p > 1)) 
    stop(paste("p must be between 0 and 1", "\n", ""))
  
  q <- ((-1/mu)*log(p))^(-1/sigma)
  q
}
#' @importFrom stats runif
#' @export
#' @rdname dIW
rIW <- function(n,mu,sigma){
  if(any(n <= 0))
    stop(paste("n must be positive","\n",""))
  if (any(mu <= 0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  
  n <- ceiling(n)
  p <- runif(n)
  r <- qIW(p, mu,sigma)
  r
}
#' @export
#' @rdname dIW
hIW<-function(x, mu, sigma){
  if (any(x < 0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(mu <= 0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  
  h <- dIW(x, mu, sigma, log=FALSE) / 
    pIW(q=x, mu, sigma, lower.tail=FALSE, log.p=FALSE)
  h
}
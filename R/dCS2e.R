#' The Cosine Sine Exponential distribution
#' 
#' @author Juan Pablo Ramirez
#' 
#' @description 
#' Density, distribution function, quantile function, 
#' random generation and hazard function for the Cosine Sine Exponential distribution 
#' with parameters \code{mu}, \code{sigma} and \code{nu}.
#' 
#' @param x,q	vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. 
#' @param mu parameter.
#' @param sigma parameter.
#' @param nu parameter.
#' @param log,log.p	logical; if TRUE, probabilities p are given as log(p).	
#' @param lower.tail logical; if TRUE (default), probabilities are 
#' P[X <= x], otherwise, P[X > x].
#' 
#' @seealso \link{CS2e}
#' 
#' @details 
#' The Cosine Sine Exponential Distribution with parameters \code{mu}, 
#' \code{sigma} and \code{nu} has density given by
#' 
#' \eqn{f(x)=\frac{\pi \sigma \mu \exp(\frac{-x} {\nu})}{2 \nu [(\mu\sin(\frac{\pi}{2} \exp(\frac{-x} {\nu})) + \sigma\cos(\frac{\pi}{2} \exp(\frac{-x} {\nu}))]^2}, }
#' 
#' for \eqn{x > 0}, \eqn{\mu > 0}, \eqn{\sigma > 0} and \eqn{\nu > 0}. 
#' 
#' @return 
#' \code{dCS2e} gives the density, \code{pCS2e} gives the distribution 
#' function, \code{qCS2e} gives the quantile function, \code{rCS2e}
#' generates random deviates and \code{hCS2e} gives the hazard function.
#'
#' @example examples/examples_dCS2e.R
#'       
#' @references
#' Chesneau, C., Bakouch, H. S., & Hussain, T. (2019). A new class 
#' of probability distributions via cosine and sine functions 
#' with applications. Communications in Statistics-Simulation 
#' and Computation, 48(8), 2287-2300.
#'
#' @export
dCS2e <- function(x, mu, sigma, nu, log=FALSE){
  if (any(x < 0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(mu <= 0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu < 0)) 
    stop(paste("nu must be positive", "\n", ""))
  
  term <- exp(-(x/nu))
  loglik <- log( pi* sigma * mu *term)-log((2*nu) * (mu*sin((pi/2) * term) + sigma*cos((pi/2) * term))^2)
  
  if (log == FALSE) 
    density <- exp(loglik)
  else 
    density <- loglik
  return(density) 
}
#' @export
#' @rdname dCS2e
pCS2e <- function(q, mu, sigma, nu, 
                  lower.tail=TRUE, log.p=FALSE){
  if (any(q < 0)) 
    stop(paste("q must be positive", "\n", ""))
  if (any(mu < 0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu < 0)) 
    stop(paste("nu must be positive", "\n", ""))
  
  term <- exp(-(q/nu))
  cdf <- (sigma*cos((pi/2) * term))/((mu*sin((pi/2) * term))+(sigma*cos((pi/2) * term)))
  
  if (lower.tail == TRUE) cdf <- cdf
  else cdf <- 1 - cdf 
  if (log.p == FALSE) cdf <- cdf
  else cdf <- log(cdf)
  cdf
}
#' @export
#' @rdname dCS2e
qCS2e <- function(p, mu, sigma, nu, 
                  lower.tail=TRUE, log.p=FALSE){
  if (any(mu < 0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu < 0)) 
    stop(paste("nu must be positive", "\n", ""))
  
  if (log.p == TRUE) p <- exp(p)
  else p <- p
  if (lower.tail == TRUE) p <- p
  else p <- 1 - p
  if (any(p < 0) | any(p > 1)) 
    stop(paste("p must be between 0 and 1", "\n", ""))
  
  q<- -(nu)*log((2/pi)*acos((p*mu)/sqrt(p^2*mu^2 + p^2*sigma^2 - 2*p*sigma^2 + sigma^2)))
  q
}
#' @importFrom stats runif
#' @export
#' @rdname dCS2e
rCS2e<- function(n, mu, sigma,nu){
  if(any(n <= 0))
    stop(paste("n must be positive","\n",""))
  if (any(mu < 0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu < 0)) 
    stop(paste("nu must be positive", "\n", ""))
  
  n <- ceiling(n)
  p <- runif(n)
  r <- qCS2e(p, mu, sigma,nu)
  r
}
#' @export
#' @rdname dCS2e
hCS2e <- function(x, mu, sigma, nu){
  if (any(x < 0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(mu < 0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu < 0)) 
    stop(paste("nu must be positive", "\n", ""))
  
  h <- dCS2e(x, mu, sigma, nu, log=FALSE) / 
    pCS2e(x, mu, sigma, nu , lower.tail=FALSE, log.p=FALSE)
  h
}


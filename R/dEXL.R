#' The exponentiated XLindley distribution
#' 
#' @author Manuel Gutierrez Tangarife, \email{mgutierrezta@unal.edu.co}
#' 
#' @description
#' Density, distribution function, quantile function, 
#' random generation and hazard function for the exponentiated XLindley distribution with
#' parameters \code{mu} and \code{sigma}.
#' 
#' @param x,q	vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. 
#' @param mu parameter.    
#' @param sigma parameter.
#' @param log,log.p	logical; if TRUE, probabilities p are given as log(p).	
#' @param lower.tail logical; if TRUE (default), probabilities are 
#' P[X <= x], otherwise, P[X > x].
#' 
#' @references
#' Alomair, A. M., Ahmed, M., Tariq, S., Ahsan-ul-Haq, M., & Talib, J. 
#' (2024). An exponentiated XLindley distribution with properties, 
#' inference and applications. Heliyon, 10(3).
#' 
#' @seealso \link{EXL}.
#' 
#' @details 
#' The exponentiated XLindley with parameters \code{mu} and \code{sigma}
#' has density given by
#' 
#' \eqn{
#'  f(x) = \frac{\sigma\mu^2(2+\mu + x)\exp(-\mu x)}{(1+\mu)^2}\left[1-
#'  \left(1+\frac{\mu x}{(1 + \mu)^2}\right) \exp(-\mu x)\right] ^ {\sigma-1} 
#' }
#' 
#' for \eqn{x \geq 0}, \eqn{\mu \geq 0} and \eqn{\sigma \geq 0}.
#' 
#' Note: In this implementation we changed the original parameters \eqn{\delta} for \eqn{\mu}
#' and \eqn{\alpha} for \eqn{\sigma}, we did it to implement this distribution
#' within gamlss framework.
#' 
#' @return 
#' \code{dEXL} gives the density, \code{pEXL} gives the distribution 
#' function, \code{qEXL} gives the quantile function, \code{rEXL}
#' generates random deviates and \code{hEXL} gives the hazard function.
#' 
#' @example examples/examples_dEXL.R
#' 
#' @export
dEXL <- function(x, mu, sigma, log=FALSE){
  if(any(x<0))     stop("Parameter x has to be positive or zero")
  if(any(mu<0))    stop("Parameter mu has to be positive or zero")
  if(any(sigma<0)) stop("Parameter sigma has to be positive or zero")

  p1 <- (1 + ((mu*x)/(1+mu)^2))
  p2 <- exp(-mu*x)
  p3 <- (1-(p1*p2))^(sigma-1) 
  p4 <- ((sigma*(mu^2))*(2+mu+x)) / ((1+mu)^2)
  p5 <- exp(-mu*x)

  pdf <- p4*p5*p3

  if(log)
    pdf <- log(pdf)
  else
    pdf <- pdf

  return(pdf)
}
#' @export
#' @rdname dEXL
pEXL <- function(q, mu, sigma, log.p=FALSE, lower.tail=TRUE){
  if(any(q < 0))   stop(paste("q must be positive", "\n", ""))
  if(any(mu<0))    stop("Parameter mu has to be positive or zero")
  if(any(sigma<0)) stop("Parameter sigma has to be positive or zero")
  
  p1 <- (1 + ((mu * q) / (1 + mu)^2))
  p2 <- exp(-mu * q)
  p3 <- 1 - (p1 * p2)
  cdf <- p3^sigma
  
  if (lower.tail == TRUE){
    cdf <- cdf
  } else {
    cdf <- 1 - cdf
  }
  
  if (log.p == FALSE){
    cdf <- cdf
  } else {
    cdf <- log(cdf)
  }
  
  return(cdf)
}
#' @importFrom lamW lambertWm1
#' @export
#' @rdname dEXL
qEXL <- function(p, mu, sigma, lower.tail=TRUE, log.p=FALSE){
  if(any(p < 0 | p > 1)) stop("Parameter p has to be beetwen 0 and 1")
  if(any(mu<=0))         stop("Parameter mu has to be positive")
  if(any(sigma<=0))      stop("Parameter sigma has to be positive")
  p1 <- -(1+mu)^2 / mu
  p2 <- 1/mu
  temp <- (1+mu)^2  * (p^(1/sigma)-1) / exp((1+mu)^2)
  p3 <- lambertWm1(temp)
  result <- p1 - p2 * p3
  return(result)
}
#' @export
#' @rdname dEXL
rEXL <- function(n, mu, sigma){
  if(any(mu<=0))    stop("Parameter mu has to be positive ")
  if(any(sigma<=0)) stop("Parameter sigma has to be positive")
  u <- runif(n)
  return(qEXL(p=u, mu=mu, sigma=sigma))
}
#' @export
#' @rdname dEXL
hEXL<- function(x, mu, sigma, log=FALSE){
  if(any(mu<=0))    stop("Parameter mu has to be positive ")
  if(any(sigma<=0)) stop("Parameter sigma has to be positive")
  p1 <- dEXL(x,mu,sigma,log=log)
  p2 <- 1 - pEXL(x,mu,sigma, log.p=log)
  return(p1/p2)
}

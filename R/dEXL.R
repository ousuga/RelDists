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
  if(any(mu<0))    stop("Parameter mu has to be positive or zero")
  if(any(sigma<0)) stop("Parameter sigma has to be positive or zero")
  
  # Ensure same length vector
  ly    <- max(length(x), length(mu), length(sigma))
  xx    <- rep(x, length=ly)
  mu    <- rep(mu, length=ly)
  sigma <- rep(sigma, length=ly)
  
  # Temporal change for invalid x's
  xx[x < 0] <- 0.5
  xx[is.infinite(x)] <- 0.5
  
  # pdf in log-scale
  p1 <- log(sigma) + 2*log(mu) + log(2+mu+xx)
  p2 <- -mu*xx - 2*log(1+mu)
  p3 <- (sigma-1)*log(1-(1+mu*xx/(1+mu)^2)*exp(-mu*xx))
  p <- p1 + p2 + p3
  
  # Assign values for invalid x's
  p[x < 0] <- -Inf
  p[is.infinite(x)] <- -Inf
  
  if (log == FALSE)
    p <- exp(p)
  
  return(p)
}
#' @export
#' @rdname dEXL
pEXL <- function(q, mu, sigma, log.p=FALSE, lower.tail=TRUE){
  if(any(mu<0))    stop("Parameter mu has to be positive or zero")
  if(any(sigma<0)) stop("Parameter sigma has to be positive or zero")
  
  # Ensure same length vector
  ly    <- max(length(q), length(mu), length(sigma))
  qq    <- rep(q, length=ly)
  mu    <- rep(mu, length=ly)
  sigma <- rep(sigma, length=ly)
  
  # Temporal change for invalid x's
  qq[q < 0] <- 0.5
  qq[q == Inf] <- 0.5
  
  # The cumulative
  p1 <- (1 + ((mu * qq) / (1 + mu)^2))
  p2 <- exp(-mu * qq)
  p3 <- 1 - (p1 * p2)
  cdf <- p3^sigma
  
  # Assign values for invalid x's
  cdf[q < 0] <- 0
  cdf[q == Inf] <- 1
  
  if (lower.tail == FALSE)
    cdf <- 1 - cdf
  if (log.p == TRUE)
    cdf <- log(cdf)
  
  return(cdf)
}
#' @importFrom lamW lambertWm1
#' @export
#' @rdname dEXL
qEXL <- function(p, mu, sigma, lower.tail=TRUE, log.p=FALSE){
  if(any(mu<=0))         stop("Parameter mu has to be positive")
  if(any(sigma<=0))      stop("Parameter sigma has to be positive")
  
  # To adjust the probability
  if (log.p == TRUE)
    p <- exp(p)
  if (lower.tail == FALSE)
    p <- 1 - p
  
  # Ensure same length vector
  ly <- max(length(p), length(mu), length(sigma))
  pp <- rep(p, length=ly)
  mu <- rep(mu, length=ly)
  sigma <- rep(sigma, length=ly)
  
  # Temporal change for invalid p's
  pp[p < 0]  <-  0.5
  pp[p > 1]  <-  0.5
  pp[p == 1] <-  0.5
  pp[p == 0] <-  0.5
  
  # The quantile
  p1 <- -(1+mu)^2 / mu
  p2 <- 1/mu
  temp <- (1+mu)^2  * (pp^(1/sigma)-1) / exp((1+mu)^2)
  p3 <- lambertWm1(temp)
  q <- p1 - p2 * p3
  
  # To deal with invalid p's
  q[p <  0] <- NaN
  q[p >  1] <- NaN
  q[p == 1] <- Inf
  q[p == 0] <- 0
  
  return(q)
}
#' @export
#' @rdname dEXL
rEXL <- function(n, mu, sigma){
  if(any(mu<=0))    stop("Parameter mu has to be positive ")
  if(any(sigma<=0)) stop("Parameter sigma has to be positive")
  if (any(n <= 0))  stop(paste("n must be a positive integer", "\n", ""))
  
  n <- ceiling(n)
  u <- runif(n=n)
  x <- qEXL(p=u, mu=mu, sigma=sigma)
  return(x)
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

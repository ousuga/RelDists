#' The Odd Weibull Distribution
#' 
#' @description 
#' Density, distribution function, quantile function, 
#' random generation and hazard function for the odd Weibull distribution with
#' parameters \code{mu}, \code{sigma} and \code{nu}.
#' 
#' @param x,q	vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. 
#' @param mu scale parameter.
#' @param sigma,nu shape parameters.
#' @param log,log.p	logical; if TRUE, probabilities p are given as log(p).	
#' @param lower.tail logical; if TRUE (default), probabilities are P[X <= x], otherwise, P[X > x].
#' 
#' @seealso \link{OW}
#' 
#' @details 
#' The odd Weibull Distribution with parameters \code{mu}, 
#' \code{sigma} and \code{nu} has density given by
#' 
#' \eqn{f(x)=\mu \sigma \nu x^{\sigma-1} exp^{\mu x^{\sigma}} (exp^{\mu x^{\sigma}}-1)^{\nu -1}(1+(exp^{\mu x^{\sigma}}-1)^{\nu})^{-2},}
#' 
#' for x > 0.
#' 
#' @return 
#' \code{dOW} gives the density, \code{pOW} gives the distribution 
#' function, \code{qOW} gives the quantile function, \code{rOW}
#' generates random deviates and \code{hOW} gives the hazard function.
#' 
#' @examples  
#' ## The probability density function 
#' curve(dOW(x, mu = 2, sigma = 3, nu = 0.2), from = 0, to = 4, ylim = c(0, 2), col = "red", las = 1, ylab = "The probability density function")
#' 
#' ## The cumulative distribution and the Reliability function
#' par(mfrow = c(1, 2))
#' curve(pOW(x, mu = 2, sigma = 3, nu = 0.2), from = 0, to = 4, ylim = c(0, 1), col = "red", las = 1, ylab = "The cumulative distribution function")
#' curve(pOW(x, mu = 2, sigma = 3, nu = 0.2, lower.tail = FALSE), from = 0, to = 4,  ylim = c(0, 1), col = "red", las = 1, ylab = "The Reliability function")
#' 
#' ## The quantile function
#' p <- seq(from = 0, to = 0.998, length.out = 100)
#' plot(x = qOW(p, mu = 2, sigma = 3, nu = 0.2), y = p, xlab = "Quantile", las = 1, ylab = "Probability")
#' curve(pOW(x, mu = 2, sigma = 3, nu = 0.2), from = 0, add = TRUE, col = "red")
#' 
#' ## The random function
#' hist(rOW(n = 10000, mu = 2, sigma = 3, nu = 0.2), freq = FALSE, ylim = c(0, 2),xlab = "x", las = 1, main = "")
#' curve(dOW(x, mu = 2, sigma = 3, nu = 0.2),  from = 0, ylim = c(0, 2), add = TRUE, col = "red")
#' 
#' ## The Hazard function
#' curve(hOW(x, mu = 2, sigma = 3, nu = 0.2), from = 0, to = 2.5, ylim = c(0, 3), col = "red", ylab = "The hazard function", las = 1)
#'
#' @export
dOW<-function(x, mu, sigma, nu, log = FALSE){
  if (any(x<0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(mu<=0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma*nu<=0)) 
    stop(paste("Product sigma*nu must be positive", "\n", ""))
  
  loglik<- log(mu) +log(sigma*nu) + (sigma-1)*log(x) +
    mu*(x^sigma) + (nu-1)*log(exp(mu*(x^sigma))-1) -
    2*log(1+(exp(mu*(x^sigma))-1)^nu)
  
  if (log == FALSE) 
    density<- exp(loglik)
  else 
    density <- loglik
  return(density)
}

#' @export
#' @rdname dOW
pOW <- function(q,mu,sigma,nu, lower.tail=TRUE, log.p = FALSE){
  if (any(q<0)) 
    stop(paste("q must be positive", "\n", ""))
  if (any(mu<=0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma*nu<=0)) 
    stop(paste("Product sigma*nu must be positive", "\n", ""))
  
  cdf <- 1 - (1 + (exp(mu*(q^sigma))-1)^nu )^(-1)
  if (lower.tail == TRUE) 
    cdf <- cdf
  else cdf <- 1 - cdf
  if (log.p == FALSE) 
    cdf <- cdf
  else cdf <- log(cdf)
  cdf
}
#' @export
#' @rdname dOW
qOW <- function(p, mu, sigma, nu, lower.tail = TRUE, log.p = FALSE){
  if (any(mu<=0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma*nu<=0)) 
    stop(paste("Product sigma*nu must be positive", "\n", ""))
  
  if (log.p == TRUE) 
    p <- exp(p)
  else p <- p
  if (lower.tail == TRUE) 
    p <- p
  else p <- 1 - p
  if (any(p < 0) | any(p > 1)) 
    stop(paste("p must be between 0 and 1", "\n", ""))
  
  q <- (1/mu)*(log( 1 + (-1+(1-p)^(-1))^(1/nu) ))^(1/sigma)
  q
}
#' @importFrom stats runif
#' @export
#' @rdname dOW
rOW <- function(n, mu, sigma, nu){
  if(any(n<=0))
    stop(paste("n must be positive","\n",""))
  if (any(mu<=0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma*nu<=0)) 
    stop(paste("Product sigma*nu must be positive", "\n", ""))
  
  n <- ceiling(n)
  p <- runif(n)
  r <- qOW(p, mu, sigma, nu)
  r
}
#' @export
#' @rdname dOW
hOW<-function(x,mu,sigma,nu){
  if (any(x<0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(mu<=0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma*nu<=0)) 
    stop(paste("Product sigma*nu must be positive", "\n", ""))
  
  h <- dOW(x, mu, sigma, nu, log = FALSE)/pOW(q = x, mu, sigma, nu, lower.tail=FALSE, log.p = FALSE)
  h
}
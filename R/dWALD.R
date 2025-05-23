#' The Wald distribution
#'
#' @author Sofia Cuartas, \email{scuartasg@unal.edu.co}
#'
#' @description
#' These functions define the density, distribution function, quantile
#' function and random generation for the Wald distribution
#' with parameter \eqn{\mu} and \eqn{\sigma}.
#'
#' @param x,q vector of (non-negative integer) quantiles.
#' @param p vector of probabilities.
#' @param mu vector of the mu parameter.
#' @param sigma vector of the sigma parameter.
#' @param n number of random values to return.
#' @param log,log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE (default), probabilities are 
#' P[X <= x], otherwise, P[X > x].
#' 
#' @seealso \link{WALD}
#'
#' @references
#' Heathcote, A. (2004). Fitting Wald and ex-Wald distributions to 
#' response time data: An example using functions for the S-PLUS package. 
#' Behavior Research Methods, Instruments, & Computers, 36, 678-694.
#'
#' @seealso \link{WALD}.
#'
#' @details
#' The Wald distribution with parameters \eqn{\mu} and \eqn{\sigma} has density given by
#'
#' \eqn{f(x |\mu, \sigma)=\frac{\sigma}{\sqrt{2 \pi x^3}} \exp \left[-\frac{(\sigma-\mu x)^2}{2x}\right ],}
#' 
#' for \eqn{x < 0}.
#'
#' @return
#' \code{dWALD} gives the density, \code{pWALD} gives the distribution
#' function, \code{qWALD} gives the quantile function, \code{rWALD}
#' generates random deviates.
#'
#' @example  examples/examples_dWALD.R
#'
#' @export
dWALD <- function(x, mu, sigma, log=FALSE) {
  if (any(x <= -1e-15)) stop(paste("x must be positive", "\n", ""))
  if (any(mu <= 0))     stop(paste("mu must  be positive 0", "\n", ""))
  if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", ""))
  
  res <- log(sigma) - 0.5*log(2*pi) - 1.5*log(x) - (sigma-mu*x)^2/(2*x)
  
  if(log == FALSE)
    return(exp(res))
  else
    return(res)
}
#' @export
#' @rdname dWALD
#' @importFrom stats pnorm
pWALD <- function(q, mu , sigma, lower.tail = TRUE, log.p = FALSE){
  if (any(q <= -1e-15)) stop(paste("q must be positive", "\n", ""))
  if (any(mu <= 0))     stop(paste("mu must  be positive", "\n", ""))
  if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", ""))
  sqrtq <- sqrt(q)
  k1 <- (mu*q-sigma)/sqrtq
  k2 <- (mu*q+sigma)/sqrtq
  p1 <- exp(2*sigma*mu)
  p2 <- pnorm(-k2)
  bad <- (p1==Inf) | (p2==0)
  p <- p1*p2
  p[bad] <- (exp(-(k1[bad]^2)/2 - 0.94/(k2[bad]^2))/(k2[bad]*((2*pi)^.5)))
  cdf <- p + pnorm(k1)
  if (lower.tail == TRUE) {
    cdf <- cdf
  }
  else {
    cdf <- 1 - cdf
  }
  if (log.p == FALSE){
    cdf <- cdf}
  else {cdf <- log(cdf)}
  return(cdf)
}
#' @export
#' @rdname dWALD
qWALD <- function(p, mu, sigma, lower.tail=TRUE, log.p=FALSE){
  if (any(mu <= 0))    stop(paste("mu must  be positive 0", "\n", ""))
  if (any(sigma <= 0)) stop(paste("sigma must be positive", "\n", ""))
  
  if (log.p == TRUE) p <- exp(p)
  else p <- p
  if (lower.tail == TRUE) p <- p
  else p <- 1 - p
  
  if (any(p < 0) | any(p > 1))  stop(paste("p must be between 0 and 1", "\n", ""))
  
  F.inv <- function(y, mu, sigma, nu) {
    uniroot(function(x) {pWALD(x,mu,sigma) - y},
            interval=c(0, 99999))$root
  }
  F.inv <- Vectorize(F.inv)
  F.inv(p, mu, sigma)
}
#' @export
#' @rdname dWALD
#' @importFrom stats rchisq
rWALD <- function(n, mu, sigma){
  if (any(n <= 0))     stop(paste("n must be positive","\n",""))
  if (any(mu <= 0))    stop(paste("mu must  be positive 0", "\n", ""))
  if (any(sigma <= 0)) stop(paste("sigma must be positive", "\n", ""))
  y2 <- rchisq(n, 1)
  y2onm <- y2/mu
  u <- runif(n)
  r1 <- (2*sigma + y2onm - sqrt(y2onm*(4*sigma+y2onm)))/(2*mu)
  r2 <- (sigma/mu)^2/r1
  r <- ifelse(u < sigma/(sigma+mu*r1), r1, r2)
  r
}

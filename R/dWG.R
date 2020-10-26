#' The Weibull Geometric distribution
#' 
#' @author Johan David Marin Benjumea, \email{johand.marin@@udea.edu.co}
#' 
#' @description 
#' Density, distribution function, quantile function, 
#' random generation and hazard function for the weibull geometric distribution with
#' parameters \code{mu}, \code{sigma} and \code{nu}.
#' 
#' @param x,q	vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. 
#' @param mu scale parameter.
#' @param sigma shape parameter.
#' @param nu parameter of geometric random variable.             
#' @param log,log.p	logical; if TRUE, probabilities p are given as log(p).	
#' @param lower.tail logical; if TRUE (default), probabilities are 
#' P[X <= x], otherwise, P[X > x].
#'  
#' @details 
#' The Weibull geometric distribution with parameters \code{mu},
#' \code{sigma} and \code{nu} has density given by
#' 
#' \eqn{f(x) = (\sigma \mu^\sigma (1-\nu) x^(\sigma - 1) \exp(-(\mu x)^\sigma)) 
#' (1- \nu \exp(-(\mu x)^\sigma))^{-2},}
#' 
#' for \eqn{x > 0}, \eqn{\mu > 0}, \eqn{\sigma > 0} and \eqn{0 < \nu < 1}.
#'
#' @return 
#' \code{dWG} gives the density, \code{pWG} gives the distribution 
#' function, \code{qWG} gives the quantile function, \code{rWG}
#' generates random deviates and \code{hWG} gives the hazard function.
#' 
#' @example examples/examples_dWG.R
#' 
#' @references
#'\insertRef{barreto2011weibull}{RelDists}
#'
#' @export
dWG<-function(x, mu, sigma, nu, log=FALSE){
  if (any(x < 0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(sigma <= 0 )) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(mu <= 0)) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(nu <= 0  | nu >= 1  )) 
    stop(paste("nu must be between zero and one", "\n", ""))
  
  loglik <- log(sigma) + sigma*log(mu) + log(1-nu) + (sigma-1)*log(x) - 
    (mu*x)^sigma - 2*log(1-nu*exp(-(mu*x)^sigma))
  
  if (log == FALSE) density<- exp(loglik)
  else density <- loglik
  return(density)
}
#' @export
#' @rdname dWG
pWG <- function(q, mu, sigma, nu, lower.tail=TRUE, log.p=FALSE){
  if (any(q < 0)) 
    stop(paste("q must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(mu <= 0)) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(nu <= 0  | nu >= 1  )) 
    stop(paste("p must be between zero and one", "\n", ""))
  
  cdf <- (1 - exp(-(mu*q)^sigma)) / (1 - nu*exp(-(mu*q)^sigma))
  
  if (lower.tail == TRUE) cdf <- cdf
  else cdf <- 1 - cdf
  if (log.p == FALSE) cdf <- cdf
  else cdf <- log(cdf)
  cdf
}
#' @export
#' @rdname dWG
qWG <- function(p, sigma, mu, nu, lower.tail = TRUE, log.p = FALSE) {
  if (any(sigma <= 0 )) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(mu <= 0)) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(nu <= 0  | nu >= 1  )) 
    stop(paste("nu must be between zero and one", "\n", ""))
  
  if (log.p == TRUE) p <- exp(p)
  else p <- p
  if (lower.tail == TRUE) p <- p
  else p <- 1 - p
  if (any(p < 0) | any(p > 1)) 
    stop(paste("p must be between 0 and 1", "\n", ""))
  
  fda <- function(x,sigma, mu,nu){
    (1- exp(-(mu*x)^sigma))/(1-(nu*exp(-(mu*x)^sigma)))
  }
  fda1 <- function(x, sigma, mu, nu, p) {fda(x, sigma, mu,nu) - p}
  r_de_la_funcion <- function(sigma, mu, nu,p) {
    uniroot(fda1, interval=c(0,1e+06), sigma, mu, nu,p)$root
  }
  r_de_la_funcion <- Vectorize(r_de_la_funcion)
  q <- r_de_la_funcion(sigma, mu, nu,p)
  q
}
#' @importFrom stats runif
#' @export
#' @rdname dWG
rWG <- function(n, mu, sigma, nu){
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(mu <= 0)) 
    stop(paste("mu must be positive", "\n", ""))  
  if (any(nu <= 0  | nu >= 1)) 
    stop(paste("nu must be between zero and one", "\n", ""))
  
  n <- ceiling(n)
  p <- runif(n)
  r <- qWG(p, mu, sigma, nu)
  r
}
#' @export
#' @rdname dWG
hWG<-function(x, mu, sigma, nu){
  if (any(x < 0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(sigma <= 0 )) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(mu <= 0)) 
    stop(paste("mu must be positive", "\n", ""))  
  if (any(nu <= 0  | nu >= 1)) 
    stop(paste("nu must be between zero and one", "\n", ""))
  
  h <- dWG(x,mu, sigma,nu, log = FALSE) / 
    pWG(q=x,mu, sigma,nu, lower.tail=FALSE, log.p = FALSE)
  h  
}

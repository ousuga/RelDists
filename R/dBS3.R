#' The Birnbaum-Saunders distribution - Bourguignon & Gallardo (2022)
#' 
#' @description
#' Density, distribution function, quantile function, 
#' random generation and hazard function for the 
#' Birnbaum-Saunders distribution with
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
#' Bourguignon, M., & Gallardo, D. I. (2022). A new look at the 
#' Birnbaum–Saunders regression model. Applied Stochastic Models in 
#' Business and Industry, 38(6), 935-951.
#' 
#' @seealso \link{BS3}.
#' 
#' @details 
#' The Birnbaum-Saunders with parameters \code{mu} and \code{sigma}
#' has density given by
#' 
#' \eqn{f(x) = \frac{\exp(\sigma/2)\sqrt{\sigma+1}}{4\sqrt{\pi\mu}x^{3/2}} \left[ x + \frac{\mu\sigma}{\sigma+1} \right] \exp\left( \frac{-\sigma}{4} \left(\frac{x(\sigma+1)}{\mu\sigma}+\frac{\mu\sigma}{x(\sigma+1)} \right) \right) }
#' 
#' for \eqn{x>0}, \eqn{\mu>0} and \eqn{\sigma>0}. In this 
#' parameterization 
#' \eqn{E(X)=\mu} and 
#' \eqn{Var(X)=(\mu\sigma)^2(1+5\sigma^2/4)}. The functions 
#' proposed here corresponds to the 
#' parameterization proposed by 
#' Santos-Neto et al. (2014).
#' 
#' @return 
#' \code{dBS3} gives the density, \code{pBS3} gives the distribution 
#' function, \code{qBS3} gives the quantile function, \code{rBS3}
#' generates random deviates and \code{hBS3} gives the hazard function.
#' 
#' @example examples/examples_dBS3.R
#' 
#' @export
dBS3 <- function(x, mu=1, sigma=0.5, log=FALSE){
  if (any(mu <= 0)) stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0 | sigma>=1)) 
    stop(paste("sigma must be in (0, 1)", "\n", ""))
  
  # Changing from BS3 to BS (original)
  new_mu    <- mu / (1-sigma)
  new_sigma <- sqrt(sigma)
  
  res <- dBS(x=x, mu=new_mu, sigma=new_sigma, log=log)
  return(res)
}
#' @export
#' @importFrom stats pnorm
#' @rdname dBS3
pBS3 <- function(q, mu=1, sigma=0.5, lower.tail=TRUE, log.p=FALSE){
  if (any(mu <= 0))    stop("parameter mu has to be positive!")
  if (any(sigma <= 0 | sigma>=1)) 
    stop(paste("sigma must be in (0, 1)", "\n", ""))
  
  # Changing from BS3 to BS (original)
  new_mu    <- mu / (1-sigma)
  new_sigma <- sqrt(sigma)
  
  cdf <- pBS(q=q, mu=new_mu, sigma=new_sigma, lower.tail=lower.tail, log.p=log.p)

  return(cdf)
}
#' @importFrom stats uniroot qnorm
#' @export
#' @rdname dBS3
qBS3 <- function(p, mu=1, sigma=0.5, lower.tail = TRUE, log.p = FALSE){
  if (any(mu <= 0)) stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0 | sigma>=1)) 
    stop(paste("sigma must be in (0, 1)", "\n", ""))
  
  # Changing from BS3 to BS (original)
  new_mu    <- mu / (1-sigma)
  new_sigma <- sqrt(sigma)
  
  if (log.p==TRUE) p <- log(p)
  if (lower.tail==FALSE) p <- 1-p
  if (any(p < 0)|any(p > 1)) stop(paste("p must be between 0 and 1", "\n", ""))

  q <- qBS(p=p, mu=new_mu, sigma=new_sigma, lower.tail=lower.tail, log.p=log.p)
  return(q)
}
#' @importFrom stats runif
#' @export
#' @rdname dBS3
rBS3 <- function(n, mu=1, sigma=0.5){
  if (any(n <= 0)) stop(paste("n must be a positive integer", "\n", ""))
  if (any(mu <= 0)) stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0 | sigma>=1)) 
    stop(paste("sigma must be in (0, 1)", "\n", ""))
  
  # Changing from BS3 to BS (original)
  new_mu    <- mu / (1-sigma)
  new_sigma <- sqrt(sigma)
  
  r <- rBS(n=n, mu=new_mu, sigma=new_sigma)
  r
}
#' @export
#' @rdname dBS3
hBS3 <- function(x, mu, sigma){
  if (any(x < 0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(mu <= 0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0 | sigma>=1)) 
    stop(paste("sigma must be in (0, 1)", "\n", ""))
  
  h <- dBS3(x, mu, sigma) / pBS3(x, mu, sigma, lower.tail=FALSE)
  h
}
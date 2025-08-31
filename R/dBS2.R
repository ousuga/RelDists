#' The Birnbaum-Saunders distribution - Santos-Neto et al. (2014)
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
#' Santos-Neto, M., Cysneiros, F. J. A., Leiva, V., & Barros, M. 
#' (2014). A reparameterized Birnbaum–Saunders distribution 
#' and its moments, estimation and applications. 
#' REVSTAT-Statistical Journal, 12(3), 247-272.
#' 
#' @seealso \link{BS2}.
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
#' \code{dBS2} gives the density, \code{pBS2} gives the distribution 
#' function, \code{qBS2} gives the quantile function, \code{rBS2}
#' generates random deviates and \code{hBS2} gives the hazard function.
#' 
#' @example examples/examples_dBS2.R
#' 
#' @export
dBS2 <- function(x, mu=1, sigma=1, log=FALSE){
  if (any(mu <= 0)) stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) stop(paste("sigma must be positive", "\n", ""))
  
  # Begin auxiliary function
  aux_fun <- function(x, mu, sigma) {
    part1 <- sigma/2 + 0.5*log(sigma+1) - log(4) - 0.5*log(pi*mu)
    part2 <- -1.5*log(x) + log(x+mu*sigma/(sigma+1))
    part3 <- (-sigma/4)*(x*(sigma+1)/(mu*sigma) + mu*sigma/(x*(sigma+1)))
    return(part1 + part2 + part3)
  }
  # End auxiliary function
  
  res <- ifelse(x<=0, -9999999, aux_fun(x, mu, sigma))
  
  if (log == TRUE)
    result <- res
  else
    result <- exp(res)
  return(result)
}
#' @export
#' @importFrom stats pnorm
#' @rdname dBS2
pBS2 <- function(q, mu=1, sigma=1, lower.tail=TRUE, log.p=FALSE){
  
  if (any(mu <= 0))    stop("parameter mu has to be positive!")
  if (any(sigma <= 0)) stop("parameter sigma has to be positive!")
  
  cdf <- pnorm(sqrt(sigma/2) * (sqrt(q*(sigma+1)/(mu*sigma)) - sqrt(mu*sigma/(q*(sigma+1)))))
  
  if (lower.tail == TRUE) 
    cdf <- cdf
  else 
    cdf = 1 - cdf
  if (log.p == FALSE) 
    cdf <- cdf
  else 
    cdf <- log(cdf)
  cdf <- ifelse(q < 0, 0, cdf)
  return(cdf)
}
#' @importFrom stats uniroot qnorm
#' @export
#' @rdname dBS2
qBS2 <- function(p, mu=1, sigma=1, lower.tail = TRUE, log.p = FALSE){
  if (any(mu <= 0)) stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) stop(paste("sigma must be positive", "\n", ""))
  if (log.p==TRUE) p <- log(p)
  if (lower.tail==FALSE) p <- 1-p
  if (any(p < 0)|any(p > 1)) stop(paste("p must be between 0 and 1", "\n", ""))
  z <- qnorm(p)   # Standard normal quantile
  term <- z / sqrt(2 * sigma)
  q <- (sigma * mu) / (sigma + 1) * (term + sqrt(term^2 + 1))^2
  return(q)
}
#' @importFrom stats runif
#' @export
#' @rdname dBS2
rBS2 <- function(n, mu=1, sigma=1){
  if (any(mu <= 0)) stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) stop(paste("sigma must be positive", "\n", ""))
  if (any(n <= 0)) stop(paste("n must be a positive integer", "\n", ""))
  n <- ceiling(n)
  p <- runif(n)
  r <- qBS2(p, mu=mu, sigma=sigma)
  r
}
#' @export
#' @rdname dBS2
hBS2 <- function(x, mu, sigma){
  if (any(x < 0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(mu <= 0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  
  h <- dBS2(x, mu, sigma) / pBS2(x, mu, sigma, lower.tail=FALSE)
  h
}
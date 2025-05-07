#' The two-parameter Chris-Jerry distribution
#' 
#' @author Manuel Gutierrez Tangarife, \email{mgutierrezta@unal.edu.co}
#' 
#' @description
#' Density, distribution function, quantile function, 
#' random generation and hazard function for the  two-parameter 
#' Chris-Jerry distribution with
#' parameters \code{mu} and \code{sigma}.
#' 
#' @param x,q	vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. 
#' @param mu parameter.    
#' @param sigma parameter.
#' @param log,log.p	logical; if TRUE, probabilities p are given as log(p).	
#' @param lower.tail logical; if TRUE (default), probabilities 
#' are P[X <= x], otherwise, P[X > x].
#' 
#' @seealso \link{CJ2}
#' 
#' @references
#' Chinedu, Eberechukwu Q., et al. "New lifetime distribution with applications 
#' to single acceptance sampling plan and scenarios of increasing hazard 
#' rates" Symmetry 15.10 (2023): 188.
#' 
#' @details 
#' The two-parameter Chris-Jerry distribution with parameters \code{mu} 
#' and \code{sigma} has density given by
#' 
#' \eqn{
#' f(x; \sigma, \mu) = \frac{\mu^2}{\sigma \mu + 2} (\sigma + \mu x^2) e^{-\mu x}; \quad x > 0, \quad \mu > 0, \quad \sigma > 0
#' }
#' 
#' Note: In this implementation we changed the original parameters 
#' \eqn{\theta} for \eqn{\mu} and \eqn{\lambda} for \eqn{\sigma},
#' we did it to implement this distribution within gamlss framework.
#' 
#' @return 
#' \code{dCJ2} gives the density, \code{pCJ2} gives the distribution 
#' function, \code{qCJ2} gives the quantile function, \code{rCJ2}
#' generates random deviates and \code{hCJ2} gives the hazard function.
#' 
#' @example examples/examples_dCJ2.R
#' 
#' @export
dCJ2 <- function(x, mu, sigma, log=FALSE){
  if(any(x<=0))     stop("Parameter x has to be positive or zero")
  if(any(mu<=0))    stop("Parameter mu has to be positive or zero")
  if(any(sigma<=0)) stop("Parameter sigma has to be positive or zero")
  
  p1 <- (mu^2/(mu*sigma + 2))
  p2 <- (sigma + mu*x^2) 
  p3 <- exp(-mu*x)
  
  pdf <- p1*p2*p3
  
  if(log)
    pdf <- log(pdf)
  else
    pdf <- pdf
  
  return(pdf)
}
#' @export
#' @rdname dCJ2
pCJ2 <- function(q, mu, sigma, log.p=FALSE, lower.tail=TRUE){
  if(any(q < 0))   stop(paste("q must be positive", "\n", ""))
  if(any(mu<=0))    stop("Parameter mu has to be positive or zero")
  if(any(sigma<=0)) stop("Parameter sigma has to be positive or zero")
  
  p1 <- (1/(mu*sigma + 2))
  p2 <- (mu^2)*(q^2) + 2*mu*q + mu*sigma + 2
  p3 <- exp(-mu * q)
  
  cdf <- 1 - (p1*p2*p3) 
  
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
#' @export
#' @rdname dCJ2
qCJ2 <- function(p, mu, sigma, lower.tail = TRUE, log.p = FALSE) {
  if (log.p) p <- exp(p)
  if (!lower.tail) p <- 1 - p
  
  if (any(p < 0 | p > 1)) stop("The probabilities 'p' must be in 0 and 1.")
  if (any(mu <= 0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  
  n <- max(length(p), length(mu), length(sigma))
  
  if (length(p)     != n) p     <- rep(p,     length.out = n)
  if (length(mu)    != n) mu    <- rep(mu,    length.out = n)
  if (length(sigma) != n) sigma <- rep(sigma, length.out = n)
  
  quant <- numeric(n)
  
  for (i in seq_len(n)) {
    root_fun <- function(x) pCJ2(x, mu[i], sigma[i]) - p[i]
    quant[i] <- uniroot(root_fun, lower = 0, upper = 1e5, tol = 1e-10)$root
  }
  
  return(quant)
}
#' @export
#' @rdname dCJ2
rCJ2 <- function(n, mu, sigma){
  if(any(mu<=0))    stop("Parameter mu has to be positive ")
  if(any(sigma<=0)) stop("Parameter sigma has to be positive")
  u <- runif(n)
  return(qCJ2(p=u, mu=mu, sigma=sigma))
}
#' @export
#' @rdname dCJ2
hCJ2 <- function(x, mu, sigma, log=FALSE){
  p1 <- dCJ2(x, mu, sigma, log=log)
  p2 <- 1 - pCJ2(x, mu, sigma, log.p=log)
  return(p1/p2)
}


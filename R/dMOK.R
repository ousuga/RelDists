#' The Marshall-Olkin Kappa distribution
#' 
#' @author Angel Muñoz,
#' 
#' @description
#' Desnsity, distribution function, quantile function,
#' random generation and hazard function for the Marshall-Olkin Kappa distribution
#' with parameters \code{mu}, \code{sigma}, \code{nu} and \code{tau}.
#' 
#' @param x,q	vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. 
#' @param mu parameter.
#' @param sigma parameter.
#' @param nu parameter.
#' @param tau parameter.
#' @param log,log.p	logical; if TRUE, probabilities p are given as log(p).	
#' @param lower.tail logical; if TRUE (default), probabilities are 
#' P[X <= x], otherwise, P[X > x].
#' 
#' @seealso \link{MOK}
#' 
#' @details
#' The Marshall-Olkin Kappa distribution with parameters \code{mu},
#' \code{sigma}, \code{nu} and \code{tau} has density given by:
#' 
#' \eqn{f(x)=\frac{\tau\frac{\mu\nu}{\sigma}\left(\frac{x}{\sigma}\right)^{\nu-1} \left(\mu+\left(\frac{x}{\sigma}\right)^{\mu\nu}\right)^{-\frac{\mu+1}{\mu}}}{\left[\tau+(1-\tau)\left(\frac{\left(\frac{x}{\sigma}\right)^{\mu\nu}}{\mu+\left(\frac{x}{\sigma}\right)^{\mu\nu}}\right)^{\frac{1}{\mu}}\right]^2}}
#' 
#' for \eqn{x > 0}.
#' 
#' @return
#' \code{dMOK} gives the density, \code{pMOK} gives the distribution function,
#' \code{qMOK} gives the quantile function, \code{rMOK} generates random deviates
#' and \code{hMOK} gives the hazard function.
#' 
#' @example examples/examples_dMOK.R
#' 
#' @references
#' Javed, M., Nawaz, T., & Irfan, M. (2019). The Marshall-Olkin 
#' kappa distribution: properties and applications. 
#' Journal of King Saud University-Science, 31(4), 684-691.
#'
#' @export
dMOK <- function(x, mu, sigma, nu, tau, log = FALSE){
  if (any(x < 0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(mu <= 0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu <= 0)) 
    stop(paste("nu must be positive", "\n", ""))
  if (any(tau <= 0)) 
    stop(paste("tau must be positive", "\n", ""))
  
  argterm <- (x/sigma)^(mu*nu)
  loglik1 <- log((tau*mu*nu)/sigma) + (nu - 1)*log(x/sigma) - ((mu + 1)/mu)*log(mu + argterm)
  loglik2 <- 2*log(tau + (1 - tau)*((argterm/(mu + argterm))^(1/mu)))
  loglik <- loglik1 - loglik2
  
  if (log == F){
    density <- exp(loglik)
  } else {
    density <- loglik
  }
  return(density)
}
#' @export
#' @rdname dMOK
pMOK <- function(q, mu, sigma, nu, tau, lower.tail = TRUE, log.p = FALSE){
  if (any(q < 0)) 
    stop(paste("q must be positive", "\n", ""))
  if (any(mu <= 0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu <= 0)) 
    stop(paste("nu must be positive", "\n", ""))
  if (any(tau <= 0)) 
    stop(paste("tau must be positive", "\n", ""))
  
  p <- (tau*((mu/((q/sigma)^(mu*nu)) + 1)^(1/mu)) + (1 - tau))^(-1)
  
  if (lower.tail == TRUE) {
    p <- p
  } else {
    p <- 1 - p
  }
  if(log.p == FALSE){
    p <- p
  } else {
    p <- log(p)
  }
  return(p)
}
#' @export
#' @rdname dMOK
qMOK <- function(p, mu, sigma, nu, tau, lower.tail = TRUE, log.p = FALSE){
  if (any(mu <= 0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu <= 0)) 
    stop(paste("nu must be positive", "\n", ""))
  if (any(tau <= 0)) 
    stop(paste("tau must be positive", "\n", ""))
  
  if (log.p == TRUE){
    p <- exp(p)
  } else {
    p <- p
  }
  if (lower.tail == TRUE){
    p <- p
  } else {
    p <- 1 - p
  }
  if (any(p < 0) | any(p > 1)) 
    stop(paste("p must be between zero and one", "\n", ""))
  
  q <- ((mu*(sigma^(mu*nu)))/(((1/p - 1 + tau)/tau)^mu - 1))^(1/(mu*nu))
  q
}
#' @export
#' @rdname dMOK
rMOK <- function(n, mu, sigma, nu, tau){
  if (any(n < 0)) 
    stop(paste("n must be positive and integer", "\n", ""))
  if (any(mu <= 0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu <= 0)) 
    stop(paste("nu must be positive", "\n", ""))
  if (any(tau <= 0)) 
    stop(paste("tau must be positive", "\n", ""))
  
  u <- runif(n)
  inv <- 1/u
  inv.a <- 1/mu
  inv.d <- 1/tau
  x <- sigma/(inv.a*((inv.d*(inv - (1 - tau)))^mu - 1))^(1/(mu*nu))
  return(x)
}
#' @export
#' @rdname dMOK
hMOK <- function(x, mu, sigma, nu, tau){
  if (any(x < 0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(mu <= 0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu <= 0)) 
    stop(paste("nu must be positive", "\n", ""))
  if (any(tau <= 0)) 
    stop(paste("tau must be positive", "\n", ""))
  
  numerator <- ((tau*mu*nu)/sigma)*((x/sigma)^(nu - 1))*((mu + (x/sigma)^(mu*nu))^(-(mu + 1)/mu))
  denominator1 <- (1 - (tau*(mu/((x/sigma)^(mu*nu)) + 1)^(1/mu) + (1 - tau))^(-1))
  denominator2 <- (tau + (1 - tau)*((((x/sigma)^(mu*nu))/(mu + (x/sigma)^(mu*nu)))^(1/mu)))^2
  h <- numerator/(denominator1*denominator2)
  return(h)
}

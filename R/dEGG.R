#' The four parameter Exponentiated Generalized Gamma distribution
#' 
#' @author Amylkar Urrea Montoya, \email{amylkar.urrea@@udea.edu.co}
#' 
#' @description 
#' Density, distribution function, quantile function, 
#' random generation and hazard function for the four parameter Exponentiated Generalized Gamma distribution
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
#' @param lower.tail logical; if TRUE (default), probabilities are P[X <= x], otherwise, P[X > x].
#' 
#' @details 
#' Four-Parameter Exponentiated Generalized Gamma distribution with parameters \code{mu}, 
#' \code{sigma}, \code{nu} and \code{tau} has density given by
#' 
#' \eqn{f(x) = \frac{\nu \sigma}{\mu \Gamma(\tau)} \left(\frac{x}{\mu}\right)^{\sigma \tau -1} \exp\left\{ - \left( \frac{x}{\mu} \right)^\sigma \right\} \left\{ \gamma_1\left( \tau, \left( \frac{x}{\mu} \right)^\sigma \right) \right\}^{\nu-1} ,}
#' 
#' for x > 0. 
#' 
#' @return 
#' \code{dEGG} gives the density, \code{pEGG} gives the distribution 
#' function, \code{qEGG} gives the quantile function, \code{rEGG}
#' generates random deviates and \code{hEGG} gives the hazard function.
#'
#' @example examples/examples_dEGG.R  
#'
#' @references
#' Almalki, S. J., & Nadarajah, S. (2014). Modifications of the 
#' Weibull distribution: A review. Reliability Engineering & 
#' System Safety, 124, 32-55.
#'
#' Cordeiro, G. M., Ortega, E. M., & Silva, G. O. (2011). 
#' The exponentiated generalized gamma distribution with 
#' application to lifetime data. Journal of statistical 
#' computation and simulation, 81(7), 827-842.
#'
#' @export
dEGG <- function(x, mu, sigma,
                   nu, tau, log=FALSE){
  if (any(x <= 0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(mu <= 0)) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", "")) 
  if (any(nu <= 0)) 
    stop(paste("nu must be positive", "\n", "")) 
  if (any(tau <= 0)) 
    stop(paste("tau must be postive", "\n", "")) 
  
  A <- log(nu * sigma) - log(mu * base::gamma(tau))
  B <- (sigma * tau - 1) * log(x / mu) - (x / mu)^sigma
  t <- (x / mu)^sigma
  C <- (nu - 1) * log(zipfR::Igamma(tau, t) / base::gamma(tau))
  loglik <- A + B + C
  
  if (log == FALSE) 
    density <- exp(loglik)
  else 
    density <- loglik
  return(density)
}
#' @export
#' @rdname dEGG
pEGG <- function(q, mu, sigma, nu, tau, 
                   lower.tail=TRUE, log.p=FALSE){
  if (any(mu <= 0)) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", "")) 
  if (any(nu <= 0)) 
    stop(paste("nu must be positive", "\n", "")) 
  if (any(tau <= 0)) 
    stop(paste("tau must be postive", "\n", "")) 
  
  
  t <- (q / mu)^sigma
  cdf <- (zipfR::Igamma(tau, t) / base::gamma(tau))^nu
  
  if (lower.tail == TRUE) 
    cdf <- cdf
  else cdf <- 1 - cdf
  if (log.p == FALSE) 
    cdf <- cdf
  else cdf <- log(cdf)
  cdf
}
#' @export
#' @rdname dEGG
qEGG <- function(p, mu, sigma, nu, tau,
                   lower.tail=TRUE, log.p=FALSE){
  if (any(mu <= 0)) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", "")) 
  if (any(nu <= 0)) 
    stop(paste("nu must be positive", "\n", "")) 
  if (any(tau <= 0)) 
    stop(paste("tau must be postive", "\n", "")) 
  if (log.p == TRUE) 
    p <- exp(p)
  else p <- p
  if (lower.tail == TRUE) 
    p <- p
  else p <- 1 - p
  if (any(p < 0) | any(p > 1)) 
    stop(paste("p must be between 0 and 1", "\n", ""))
  
  fda <- function(x, mu, sigma, nu, tau){
    
    t <- (x / mu)^sigma
    cdf <- (zipfR::Igamma(tau, t) / base::gamma(tau))^nu
    cdf
    
  }
  fda1 <- function(x, mu, sigma, nu, tau, p) {
    fda(x, mu, sigma, nu, tau) - p
  }
  r_de_la_funcion <- function(mu, sigma, nu, tau, p) {
    uniroot(fda1, interval=c(0, 1e+06), mu, sigma, nu, tau, p)$root
  }
  r_de_la_funcion <- Vectorize(r_de_la_funcion)
  q <- r_de_la_funcion(mu, sigma, nu, tau, p)
  q
}
#' @importFrom stats runif
#' @export
#' @rdname dEGG
rEGG <- function(n, mu, sigma, nu, tau){
  if (any(mu <= 0)) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", "")) 
  if (any(nu <= 0)) 
    stop(paste("nu must be positive", "\n", "")) 
  if (any(tau <= 0)) 
    stop(paste("tau must be postive", "\n", "")) 
  
  n <- ceiling(n)
  p <- runif(n)
  r <- qEGG(p, mu, sigma, nu, tau)
  r
}
#' @export
#' @rdname dEGG
hEGG <- function(x, mu, sigma, nu, tau){
  if (any(x <= 0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(mu <= 0)) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", "")) 
  if (any(nu <= 0)) 
    stop(paste("nu must be positive", "\n", "")) 
  if (any(tau <= 0)) 
    stop(paste("tau must be postive", "\n", "")) 
  
  h <- dEGG(x, mu, sigma, nu, tau, log=FALSE) / 
    pEGG(q=x, mu, sigma, nu, tau, lower.tail=FALSE, log.p=FALSE)
  h  
}

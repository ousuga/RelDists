#' The Marshall-Olkin Extended Inverse Weibull distribution
#' 
#' @author Amylkar Urrea Montoya, \email{amylkar.urrea@@udea.edu.co}
#' 
#' @description 
#' Density, distribution function, quantile function, 
#' random generation and hazard function for the Marshall-Olkin Extended Inverse Weibull distribution
#' with parameters \code{mu}, \code{sigma} and \code{nu}.
#' 
#' @param x,q	vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. 
#' @param mu parameter.
#' @param sigma parameter.
#' @param nu parameter.
#' @param log,log.p	logical; if TRUE, probabilities p are given as log(p).	
#' @param lower.tail logical; if TRUE (default), probabilities are 
#' P[X <= x], otherwise, P[X > x].
#' 
#' @seealso \link{MOEIW}
#' 
#' @details 
#' The Marshall-Olkin Extended Inverse Weibull distribution \code{mu}, 
#' \code{sigma} and \code{nu} has density given by
#' 
#' \eqn{f(x) = \frac{\mu \sigma \nu x^{-(\sigma + 1)} exp\{{-\mu x^{-\sigma}}\}}{\{\nu -(\nu-1) exp\{{-\mu x ^{-\sigma}}\} \}^{2}},}
#' 
#' for \eqn{x > 0}. 
#' 
#' @return 
#' \code{dMOEIW} gives the density, \code{pMOEIW} gives the distribution 
#' function, \code{qMOEIW} gives the quantile function, \code{rMOEIW}
#' generates random deviates and \code{hMOEIW} gives the hazard function.
#'
#' @example examples/examples_dMOEIW.R  
#'
#' @references
#' Okasha, H. M., El-Baz, A. H., Tarabia, A. M. K., & Basheer, A. M. (2017). 
#' Extended inverse Weibull distribution with reliability application. 
#' Journal of the Egyptian Mathematical Society, 25(3), 343-349.
#'
#' @export
dMOEIW <- function(x, mu, sigma, nu, log=FALSE){
  if (any(x < 0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(mu <= 0)) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", "")) 
  if (any(nu <= 0)) 
    stop(paste("nu must be positive", "\n", "")) 
  
  A <- log(mu) + log(sigma) + log(nu) - (sigma + 1) * log(x) - mu * x^(-sigma)  
  B <- 2 * log(nu - (nu - 1) * exp(-mu * x^(-sigma)))
  loglik <- A - B
  
  if (log == FALSE) 
    density <- exp(loglik)
  else 
    density <- loglik
  return(density)
}
#' @export
#' @rdname dMOEIW
pMOEIW <- function(q, mu, sigma, nu, 
                   lower.tail=TRUE, log.p=FALSE){
  if (any(mu <= 0)) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", "")) 
  if (any(nu <= 0)) 
    stop(paste("nu must be positive", "\n", "")) 
  
  A <- exp(-mu * q^(-sigma))
  B <- nu - (nu - 1) * A
  cdf <- A / B
  
  if (lower.tail == TRUE) 
    cdf <- cdf
  else cdf <- 1 - cdf
  if (log.p == FALSE) 
    cdf <- cdf
  else cdf <- log(cdf)
  cdf
}
#' @export
#' @rdname dMOEIW
qMOEIW <- function(p, mu, sigma, nu,
                   lower.tail=TRUE, log.p=FALSE){
  if (any(mu <= 0)) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", "")) 
  if (any(nu <= 0)) 
    stop(paste("nu must be positive", "\n", "")) 
  if (log.p == TRUE) 
    p <- exp(p)
  else p <- p
  if (lower.tail == TRUE) 
    p <- p
  else p <- 1 - p
  if (any(p < 0) | any(p > 1)) 
    stop(paste("p must be between 0 and 1", "\n", ""))
  
  fda <- function(x, mu, sigma, nu){
    
    exp(-mu * x^(-sigma)) / (nu - (nu - 1) * exp(-mu * x^(-sigma)))
    
  }
  fda1 <- function(x, mu, sigma, nu, p) {
    fda(x, mu, sigma, nu) - p
  }
  r_de_la_funcion <- function(mu, sigma, nu, p) {
    uniroot(fda1, interval=c(0, 1e+06), mu, sigma, nu, p)$root
  }
  r_de_la_funcion <- Vectorize(r_de_la_funcion)
  q <- r_de_la_funcion(mu, sigma, nu, p)
  q
}
#' @importFrom stats runif
#' @export
#' @rdname dMOEIW
rMOEIW <- function(n, mu, sigma, nu){
  if (any(mu <= 0)) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu <= 0)) 
    stop(paste("nu must be positive", "\n", "")) 
  
  n <- ceiling(n)
  p <- runif(n)
  r <- qMOEIW(p, mu, sigma, nu)
  r
}
#' @export
#' @rdname dMOEIW
hMOEIW<-function(x, mu, sigma, nu){
  if (any(x < 0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(mu <= 0)) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))  
  if (any(nu <= 0)) 
    stop(paste("nu must be positive", "\n", "")) 
  
  h <- dMOEIW(x, mu, sigma, nu, log=FALSE) / 
    pMOEIW(q=x, mu, sigma, nu, lower.tail=FALSE, log.p=FALSE)
  h 
}

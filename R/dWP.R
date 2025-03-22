#' The Weibull Poisson distribution
#' 
#' @author Amylkar Urrea Montoya, \email{amylkar.urrea@@udea.edu.co}
#' 
#' @description 
#' Density, distribution function, quantile function, 
#' random generation and hazard function for the Weibull Poisson distribution
#' with parameters \code{mu}, \code{sigma} and \code{nu}.
#' 
#' @param x,q	vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. 
#' @param mu parameter.
#' @param sigma parameter.
#' @param nu parameter.
#' @param log,log.p	logical; if TRUE, probabilities p are given as log(p).	
#' @param lower.tail logical; if TRUE (default), probabilities are P[X <= x], otherwise, P[X > x].
#' 
#' @details 
#' The Weibull Poisson distribution with parameters \code{mu}, 
#' \code{sigma} and \code{nu} has density given by
#' 
#' \eqn{f(x) = \frac{\mu \sigma \nu e^{-\nu}} {1-e^{-\nu}} x^{\mu-1} \exp({-\sigma x^{\mu}+\nu \exp({-\sigma} x^{\mu}) }),}
#' 
#' for \eqn{x > 0}. 
#' 
#' @return 
#' \code{dWP} gives the density, \code{pWP} gives the distribution 
#' function, \code{qWP} gives the quantile function, \code{rWP}
#' generates random deviates and \code{hWP} gives the hazard function.
#'
#' @example examples/examples_dWP.R    
#'
#' @references
#' Lu, Wanbo, and Daimin Shi. "A new compounding life distribution: 
#' the Weibull–Poisson distribution." Journal of applied 
#' statistics 39.1 (2012): 21-38.
#'
#' @importFrom Rdpack reprompt
#'
#' @export
dWP <- function(x, mu, sigma, nu, log=FALSE){
  if (any(x <= 0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(mu <= 0)) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", "")) 
  if (any(nu <= 0)) 
    stop(paste("nu must be positive", "\n", "")) 
  
  A <- log(mu) + log(sigma) + log(nu) - nu - log( 1 - exp(-nu))
  B <- (mu - 1) * log(x) + (nu * exp(-sigma * x^mu) - sigma * x^mu )
  loglik <- A + B
  
  if (log == FALSE) 
    density <- exp(loglik)
  else 
    density <- loglik
  return(density)
}
#' @export
#' @rdname dWP
pWP <- function(q, mu, sigma, nu, 
                lower.tail=TRUE, log.p=FALSE){
  if (any(mu <= 0)) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", "")) 
  if (any(nu <= 0)) 
    stop(paste("nu must be positive", "\n", "")) 
  
  A <- exp(nu * exp(- sigma * q^mu)) - exp(nu)
  cdf <- A / (1 - exp(nu))
  
  if (lower.tail == TRUE) 
    cdf <- cdf
  else cdf <- 1 - cdf
  if (log.p == FALSE) 
    cdf <- cdf
  else cdf <- log(cdf)
  cdf
}
#' @export
#' @rdname dWP
qWP <- function(p, mu, sigma, nu,
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
    
    (exp(nu * exp(- sigma * x^mu)) - exp(nu)) / (1 - exp(nu))
    
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
#' @rdname dWP
rWP <- function(n, mu, sigma, nu){
  if (any(mu <= 0)) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu <= 0)) 
    stop(paste("nu must be positive", "\n", "")) 
  
  n <- ceiling(n)
  p <- runif(n)
  r <- qWP(p, mu, sigma, nu)
  r
}
#' @export
#' @rdname dWP
hWP<-function(x, mu, sigma, nu){
  if (any(x <= 0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(mu <= 0)) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))  
  if (any(nu <= 0)) 
    stop(paste("nu must be positive", "\n", "")) 
  
  h <- dWP(x, mu, sigma, nu, log=FALSE) / 
    pWP(q=x, mu, sigma, nu, lower.tail=FALSE, log.p=FALSE)
  h  
}
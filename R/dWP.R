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
#' \eqn{f(x) = \frac{\mu \sigma \nu e^{-\nu}} {1-e^{-\nu}} x^{\mu-1} exp({-\sigma x^{\mu}+\nu exp({-\sigma} x^{\mu}) }),}
#' 
#' for x > 0. 
#' 
#' @return 
#' \code{dWP} gives the density, \code{pWP} gives the distribution 
#' function, \code{qWP} gives the quantile function, \code{rWP}
#' generates random deviates and \code{hWP} gives the hazard function.
#'
#' @examples  
#' ## The probability density function
#' curve(dWP(x, mu=1.5, sigma=0.5, nu=10), from=0.0001, to=2,
#'       col="red", las=1, ylab="f(x)")
#' 
#' ## The cumulative distribution and the Reliability function
#' par(mfrow=c(1, 2))
#' curve(pWP(x, mu=1.5, sigma=0.5, nu=10),
#'       from=0.0001, to=2, col="red", las=1, ylab="F(x)")
#' curve(pWP(x, mu=1.5, sigma=0.5, nu=10, lower.tail=FALSE),
#'       from=0.0001, to=2, col="red", las=1, ylab="S(x)")
#' 
#' ## The quantile function
#' p <- seq(from=0, to=0.99999, length.out=100)
#' plot(x=qWP(p, mu=1.5, sigma=0.5, nu=10), y=p, xlab="Quantile",
#'      las=1, ylab="Probability")
#' curve(pWP(x, mu=1.5, sigma=0.5, nu=10),
#'       from=0, add=TRUE, col="red")
#' 
#' ## The random function
#' hist(rWP(n=10000, mu=1.5, sigma=0.5, nu=10), freq=FALSE,
#'      xlab="x", ylim=c(0, 2.2), las=1, main="")
#' curve(dWP(x, mu=1.5, sigma=0.5, nu=10),
#'       from=0.001, to=4, add=TRUE, col="red")
#' 
#' ## The Hazard function
#' curve(hWP(x, mu=1.5, sigma=0.5, nu=10), from=0.001, to=5,
#'       col="red", ylab="Hazard function", las=1)
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
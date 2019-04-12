#' The Four-Parameter Generalized Gamma distribution
#' 
#' @description 
#' Density, distribution function, quantile function, 
#' random generation and hazard function for the Four-Parameter Generalized Gamma distribution
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
#' Four-Parameter Generalized Gamma distribution with parameters \code{mu}, 
#' \code{sigma}, \code{nu} and \code{tau} has density given by
#' 
#' \eqn{f(x) = \frac {\sigma \mu^{\tau \sigma}}{\tau (\tau)} (x - \nu)^{\tau \sigma -1} \exp-({\frac {(x-\nu)} {\mu} })^{\sigma},}
#' 
#' for x >= 0. 
#' 
#' @return 
#' \code{dFPGG} gives the density, \code{pFPGG} gives the distribution 
#' function, \code{qFPGG} gives the quantile function, \code{rFPGG}
#' generates random deviates and \code{hFPGG} gives the hazard function.
#'
#' @examples  
#' ## The probability density function
#' curve(dFPGG(x, mu=1, sigma=1, nu=0, tau=4), from=0.1, to=15,
#'       col="red", las=1, ylab="f(x)")
#' 
#' ## The cumulative distribution and the Reliability function
#' par(mfrow=c(1, 2))
#' curve(pFPGG(x, mu=1, sigma=1, nu=0, tau=4),
#'       from=0, to=4, col="red", las=1, ylab="F(x)")
#' curve(pFPGG(x,  mu=1, sigma=1, nu=0, tau=4, lower.tail=FALSE),
#'       from=0, to=4, col="red", las=1, ylab="S(x)")
#' 
#' ## The quantile function
#' p <- seq(from=0, to=0.99999, length.out=100)
#' plot(x=qFPGG(p, mu=1, sigma=1, nu=0, tau=4), y=p, xlab="Quantile",
#'      las=1, ylab="Probability")
#' curve(pFPGG(x,  mu=1, sigma=1, nu=0, tau=4), 
#'       from=0, add=TRUE, col="red")
#' 
## The random function
#' hist(rFPGG(n=10000,  mu=1, sigma=1, nu=0, tau=4), freq=FALSE,
#'      xlab="x", las=1, main="")
#' curve(dFPGG(x,  mu=1, sigma=1, nu=0, tau=4),
#'       from=0, to=15, add=TRUE, col="red")
#' 
#' ## The Hazard function
#' curve(hFPGG(x,  mu=1, sigma=1, nu=0, tau=4), from=0, to=30,
#'       col="red", ylab="Hazard function", las=1)
#'
#' @export
dFPGG <- function(x, mu, sigma,
                  nu, tau, log=FALSE){
  if (any(x < 0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(mu <= 0)) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", "")) 
  if (any(nu < 0)) 
    stop(paste("nu must be positive", "\n", "")) 
  if (any(tau <= 0)) 
    stop(paste("tau must be postive", "\n", "")) 
  
  A <- log(sigma) + (tau * sigma) * log(mu) - base::lgamma(tau)
  B <- (tau * sigma - 1) * log(x - nu) - ((x - nu) / mu)^sigma
  loglik <- A + B
  
  if (log == FALSE) 
    density <- exp(loglik)
  else 
    density <- loglik
  return(density)
}
#' @export
#' @rdname dFPGG
pFPGG <- function(q, mu, sigma, nu, tau, 
                  lower.tail=TRUE, log.p=FALSE){
  if (any(mu <= 0)) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", "")) 
  if (any(nu < 0)) 
    stop(paste("nu must be positive", "\n", "")) 
  if (any(tau <= 0)) 
    stop(paste("tau must be postive", "\n", "")) 
  
  
  t <- ((q - nu) / mu)^sigma
  A <- zipfR::Igamma(tau, t)
  cdf <- A / base::gamma(tau)
  
  if (lower.tail == TRUE) 
    cdf <- cdf
  else cdf <- 1 - cdf
  if (log.p == FALSE) 
    cdf <- cdf
  else cdf <- log(cdf)
  cdf
}
#' @export
#' @rdname dFPGG
qFPGG <- function(p, mu, sigma, nu, tau,
                  lower.tail=TRUE, log.p=FALSE){
  if (any(mu <= 0)) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", "")) 
  if (any(nu < 0)) 
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
    
    t <- ((x - nu) / mu)^sigma
    A <- zipfR::Igamma(tau, t)
    cdf <- A / base::gamma(tau)
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
#' @rdname dFPGG
rFPGG <- function(n, mu, sigma, nu, tau){
  if (any(mu <= 0)) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", "")) 
  if (any(nu < 0)) 
    stop(paste("nu must be positive", "\n", "")) 
  if (any(tau <= 0)) 
    stop(paste("tau must be postive", "\n", "")) 
  
  n <- ceiling(n)
  p <- runif(n)
  r <- qFPGG(p, mu, sigma, nu, tau)
  r
}
#' @export
#' @rdname dFPGG
hFPGG <- function(x, mu, sigma, nu, tau){
  if (any(x < 0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(mu <= 0)) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", "")) 
  if (any(nu < 0)) 
    stop(paste("nu must be positive", "\n", "")) 
  if (any(tau <= 0)) 
    stop(paste("tau must be postive", "\n", "")) 
  
  h <- dFPGG(x, mu, sigma, nu, tau, log=FALSE) / 
    pFPGG(q=x, mu, sigma, nu, tau, lower.tail=FALSE, log.p=FALSE)
  h  
}
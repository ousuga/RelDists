#' The Gamma Weibull distribution
#' 
#' @description 
#' Density, distribution function, quantile function, 
#' random generation and hazard function for the Gamma Weibull distribution 
#' with parameters \code{mu}, \code{sigma} and \code{nu} 
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
#' The Gamma Weibull Distribution with parameters \code{mu}, 
#' \code{sigma} and \code{nu} has density given by
#' 
#' \eqn{f(x) = (\frac{\sigma \mu ^{\nu}}{\tau (\nu)} x^{\nu \sigma - 1}) \exp({-\mu x^{\sigma} }),}
#' 
#' for x > 0. 
#' 
#' @return 
#' \code{dGW} gives the density, \code{pGW} gives the distribution 
#' function, \code{qGW} gives the quantile function, \code{rGW}
#' generates random deviates and \code{hGW} gives the hazard function.
#'
#' @examples  
#' ## The probability density function
#' curve(dGW(x, mu=1, sigma=2, nu=1), from=0.0001, to=4,
#'       col="red", las=1, ylab="f(x)")
#' 
#' ## The cumulative distribution and the Reliability function
#' par(mfrow=c(1, 2))
#' curve(pGW(x, mu=1, sigma=2, nu=1),
#'       from=0.0001, to=4, col="red", las=1, ylab="F(x)")
#' curve(pGW(x, mu=1, sigma=2, nu=1, lower.tail=FALSE),
#'       from=0.0001, to=4, col="red", las=1, ylab="S(x)")
#' 
#' ## The quantile function
#' p <- seq(from=0, to=0.99999, length.out=100)
#' plot(x=qGW(p, mu=1, sigma=2, nu=1), y=p, xlab="Quantile",
#'      las=1, ylab="Probability")
#' curve(pGW(x, mu=1, sigma=2, nu=1), 
#'       from=0.1, add=TRUE, col="red")
#' 
#' ## The random function
#' hist(rGW(n=10000, mu=1, sigma=2, nu=1), freq=FALSE,
#'      xlab="x", las=1, main="")
#' curve(dGW(x, mu=1, sigma=2, nu=1),
#'       from=0.01, to=4, add=TRUE, col="red")
#' 
#' ## The Hazard function
#' curve(hGW(x, mu=1, sigma=2, nu=1), from=1, to=6,
#'      col="red", ylab="Hazard function", las=1)
#'
#' @export
dGW <- function(x, mu, sigma, nu, log=FALSE){
  if (any(x <= 0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(mu <= 0)) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu <= 0)) 
    stop(paste("nu must be positive", "\n", "")) 
  
  loglik <- log(sigma*mu^(nu)) - log(gamma(nu)) + 
    (nu*sigma - 1)*log(x) - mu*x^(sigma)
  
  if (log == FALSE) 
    density <- exp(loglik)
  else 
    density <- loglik
  return(density)
}
#' @importFrom zipfR Igamma
#' @export
#' @rdname dGW
pGW <- function(q, mu, sigma, nu, 
                lower.tail=TRUE, log.p=FALSE){
  if (any(mu <= 0)) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", "")) 
  if (any(nu <= 0)) 
    stop(paste("nu must be positive", "\n", "")) 
  
  t <- mu * q^sigma
  A <- Igamma(nu, t)
  cdf <- A / gamma(nu)
  
  if (lower.tail == TRUE) 
    cdf <- cdf
  else cdf <- 1 - cdf
  if (log.p == FALSE) 
    cdf <- cdf
  else cdf <- log(cdf)
  cdf
}
#' @export
#' @rdname dGW
qGW <- function(p, mu, sigma, nu,
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
    t <- mu * x^sigma
    A <- Igamma(nu, t)
    cdf <- A / gamma(nu)
    cdf
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
#' @rdname dGW
rGW <- function(n, mu, sigma, nu){
  if (any(mu <= 0)) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu <= 0)) 
    stop(paste("nu must be positive", "\n", "")) 
  
  n <- ceiling(n)
  p <- runif(n)
  r <- qGW(p, mu, sigma, nu)
  r
}
#' @export
#' @rdname dGW
hGW<-function(x, mu, sigma, nu){
  if (any(x <= 0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(mu <= 0)) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))  
  if (any(nu <= 0)) 
    stop(paste("nu must be positive", "\n", "")) 
  
  h <- dGW(x, mu, sigma, nu, log=FALSE) / 
    pGW(q=x, mu, sigma, nu, lower.tail=FALSE, log.p=FALSE)
  h  
}

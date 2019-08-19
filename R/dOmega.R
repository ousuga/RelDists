#' The Omega distribution
#'
#' @author Edwin Caicedo
#'
#' @description
#' Density, distribution function, quantile function,
#' random generation and hazard function for the Omega distribution
#' with parameters \code{mu}, \code{sigma}, \code{nu}.
#'
#' @param x,q	vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations.
#' @param mu scale parameter.
#' @param sigma shape parameter.
#' @param nu parameter.
#' @param log,log.p	logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE (default), probabilities are P[X <= x], otherwise, P[X > x].
#'
#' @details
#' The Omega Distribution with parameters \code{mu},
#' \code{sigma} and \code{nu} has density given by
#'
#' \eqn{f(x)=\mu \sigma x^{\sigma - 1} \frac{\nu^{2\sigma}}{\nu^{2\sigma} - x^{2\sigma}} \omega_{\nu}^{(-\mu, \sigma)}(x),}
#'
#' where
#'
#' \eqn{\omega_{\nu}^{(-\mu, \sigma)}(x)=\left(\frac{\nu^{\sigma}+x^{\sigma}}{\nu^{\sigma}-x^{\sigma}}\right)^{\frac{-\mu \nu^{\sigma}}{2}},}
#'
#' for \eqn{x > 0}, \eqn{\mu > 0}, \eqn{\sigma > 1} and \eqn{\nu > 0}.
#'
#' @return
#' \code{dOme} gives the density, \code{pOmega} gives the distribution
#' function, \code{qOmega} gives the quantile function, \code{rOmega}
#' generates random deviates and \code{hOmega} gives the hazard function.
#' 
#' @return 
#' \code{dOmega} gives the density, \code{pOmega} gives the distribution 
#' function, \code{qOmega} gives the quantile function, \code{rOmega}
#' generates random deviates and \code{hOmega} gives the hazard function.
#' 
#' @examples
#' ## The probability density function
#' curve(dOmega(x, mu=0.4, sigma=1.5, nu=5), from=0, to=5,
#'       ylim=c(0, 0.5), col="red", las=1, ylab="f(x)")
#'
#' ## The cumulative distribution and the Reliability function
#' par(mfrow=c(1, 2))
#' curve(pOmega(x, mu=0.4, sigma=1.5, nu=5),
#'       from=0, to=2,  col="red", las=1, ylab="F(x)")
#' curve(pOmega(x, mu=0.4, sigma=1.5, nu=8, lower.tail=FALSE),
#'       from=0, to=10, col="red", las=1, ylab="S(x)")
#'
#' ## The quantile function
#' p <- seq(from=0, to=0.99999, length.out=100)
#' plot(x=qOmega(p, mu=0.3, sigma=1.5, nu=5), y=p, xlab="Quantile",
#'      las=1, ylab="Probability")
#' curve(pOmega(x, mu=0.3, sigma=1.5, nu=5), from=0, add=TRUE, col="red")
#'
#' ## The random function
#' hist(rOmega(n=10000, mu=0.4, sigma=1.5, nu=5), freq=FALSE,
#'      xlab="x", las=1, main="")
#' curve(dOmega(x, mu=0.4, sigma=1.5, nu=5), from=0, add=TRUE, col="red")
#'
#' ## The Hazard function
#' par(mfrow=c(1,1))
#' curve(hOmega(x, mu=0.4, sigma=1.5, nu=50), from=0, to=5, ylim=c(0, 2),
#'       col="red", ylab="Hazard function", las=1)
#'
#' @references
#'\insertRef{dombi2019omega}{RelDists}
#'
#' @export
dOmega <- function(x, mu, sigma, nu, log=FALSE){
  if (any(x < 0))
    stop(paste("x must be positive", "\n", ""))
  if (any(mu <= 0 ))
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0))
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu <= 0))
    stop(paste("nu must be positive", "\n", ""))
  
  term1 <- (nu^sigma + x^sigma)
  term2 <- (nu^sigma - x^sigma)
  
  loglik <- log(mu*sigma) + (sigma-1)*log(x)  + 2*sigma*log(nu) -
    log(nu^(2*sigma)-x^(2*sigma)) - (mu*(nu^sigma)/2)*log(term1/term2)
  
  if (log == FALSE)
    density <- exp(loglik)
  else
    density <- loglik
  return(density)
}
#' @export
#' @rdname dOmega
pOmega <- function(q, mu, sigma, nu,
                   lower.tail=TRUE, log.p=FALSE){
  if (any(q < 0))
    stop(paste("q must be positive", "\n", ""))
  if (any(mu <= 0 ))
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0))
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu <= 0))
    stop(paste("nu must be positive", "\n", ""))
  
  term1 <- (nu^sigma + q^sigma)
  term2 <- (nu^sigma - q^sigma)
  
  
  cdf <- 1 -  (term1/term2)^(-mu*(nu^sigma)/2)
  
  if (lower.tail == TRUE) cdf <- cdf
  else cdf <- 1 - cdf
  if (log.p == FALSE) cdf <- cdf
  else cdf <- log(cdf)
  cdf
}
#' @export
#' @rdname dOmega
qOmega <- function(p, mu, sigma, nu,
                   lower.tail=TRUE, log.p=FALSE){
  if (any(mu <= 0 ))
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0))
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu <= 0))
    stop(paste("nu must be positive", "\n", ""))
  
  if (log.p == TRUE) p <- exp(p)
  else p <- p
  if (lower.tail == TRUE) p <- p
  else p <- 1 - p
  if (any(p < 0) | any(p > 1))
    stop(paste("p must be between 0 and 1", "\n", ""))
  
  theta <- exp((-2*log(1-p))/(mu*(nu^sigma)))
  
  q <- ((theta*(nu^sigma) - nu^sigma)/(1+theta))^(1/sigma)
  q
}
#' @importFrom stats runif
#' @export
#' @rdname dOmega
rOmega <- function(n, mu, sigma, nu){
  if(any(n <= 0))
    stop(paste("n must be positive","\n",""))
  if (any(mu <= 0 ))
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0))
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu <= 0))
    stop(paste("nu must be positive", "\n", ""))
  
  n <- ceiling(n)
  p <- runif(n)
  r <- qOmega(p, mu, sigma, nu)
  r
}
#' @export
#' @rdname dOmega
hOmega <- function(x, mu, sigma, nu){
  if (any(x < 0))
    stop(paste("x must be positive", "\n", ""))
  if (any(mu <= 0 ))
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0))
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu <= 0))
    stop(paste("nu must be positive", "\n", ""))
  
  h <- dOmega(x, mu, sigma, nu, log=FALSE) /
    pOmega(x, mu, sigma, nu, lower.tail=FALSE, log.p=FALSE)
  h
}

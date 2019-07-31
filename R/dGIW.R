#' The Generalized Inverse Weibull distribution
#' 
#' @author Amylkar Urrea Montoya, \email{amylkar.urrea@@udea.edu.co}
#' 
#' @description 
#' Density, distribution function, quantile function, 
#' random generation and hazard function for the Generalized Inverse Weibull distribution
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
#' The Generalized Inverse Weibull distribution \code{mu}, 
#' \code{sigma} and \code{nu} has density given by
#' 
#' \eqn{f(x) = \nu \sigma \mu^{\sigma} x^{-(\sigma + 1)} exp \{-\nu (\frac{\mu}{x})^{\sigma}\},}
#' 
#' for x > 0. 
#' 
#' @return 
#' \code{dGIW} gives the density, \code{pGIW} gives the distribution 
#' function, \code{qGIW} gives the quantile function, \code{rGIW}
#' generates random deviates and \code{hGIW} gives the hazard function.
#'
#' @examples  
#' ## The probability density function
#' curve(dGIW(x, mu=3, sigma=5, nu=0.5), from=0.001, to=8,
#'       col="red", ylab="f(x)", las=1)
#' 
#' ## The cumulative distribution and the Reliability function
#' par(mfrow=c(1, 2))
#' curve(pGIW(x, mu=3, sigma=5, nu=0.5),
#'       from=0.0001, to=14, col="red", las=1, ylab="F(x)")
#' curve(pGIW(x, mu=3, sigma=5, nu=0.5, lower.tail=FALSE),
#'       from=0.0001, to=14, col="red", las=1, ylab="S(x)")
#' 
#' ## The quantile function
#' p <- seq(from=0, to=0.99999, length.out=100)
#' plot(x=qGIW(p, mu=3, sigma=5, nu=0.5), y=p, xlab="Quantile",
#'      las=1, ylab="Probability")
#' curve(pGIW(x, mu=3, sigma=5, nu=0.5),
#'       from=0, add=TRUE, col="red")
#' 
#' ## The random function
#' hist(rGIW(n=10000, mu=3, sigma=5, nu=0.5), freq=FALSE,
#'      xlab="x", ylim=c(0, 0.8), las=1, main="")
#' curve(dGIW(x, mu=3, sigma=5, nu=0.5),
#'       from=0.001, to=14, add=TRUE, col="red")
#' 
#' ## The Hazard function
#' curve(hGIW(x, mu=3, sigma=5, nu=0.5), from=0.001, to=30,
#'       col="red", ylab="Hazard function", las=1)
#' 
#'
#' @references
#' \insertRef{almalki2014modifications}{RelDists}
#' 
#' \insertRef{gusmao2009}{RelDists}
#'
#' @importFrom Rdpack reprompt
#'
#' @export
dGIW <- function(x, mu, sigma, nu, log=FALSE){
  if (any(x <= 0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(mu <= 0)) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", "")) 
  if (any(nu <= 0)) 
    stop(paste("nu must be positive", "\n", "")) 
  
  A <- log(nu) + log(sigma) + sigma * log(mu)
  B <- - (sigma + 1) * log(x) - nu * (mu / x)^sigma
  loglik <- A + B
  
  if (log == FALSE) 
    density <- exp(loglik)
  else 
    density <- loglik
  return(density)
}
#' @export
#' @rdname dGIW
pGIW <- function(q, mu, sigma, nu, 
                 lower.tail=TRUE, log.p=FALSE){
  if (any(mu <= 0)) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", "")) 
  if (any(nu <= 0)) 
    stop(paste("nu must be positive", "\n", "")) 
  
  cdf <- exp(-nu * (mu / q)^sigma)
  
  if (lower.tail == TRUE) 
    cdf <- cdf
  else cdf <- 1 - cdf
  if (log.p == FALSE) 
    cdf <- cdf
  else cdf <- log(cdf)
  cdf
}
#' @export
#' @rdname dGIW
qGIW <- function(p, mu, sigma, nu,
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
    
    exp(-nu * (mu / x)^sigma) 
    
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
#' @rdname dGIW
rGIW <- function(n, mu, sigma, nu){
  if (any(mu <= 0)) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu <= 0)) 
    stop(paste("nu must be positive", "\n", "")) 
  
  n <- ceiling(n)
  p <- runif(n)
  r <- qGIW(p, mu, sigma, nu)
  r
}
#' @export
#' @rdname dGIW
hGIW <- function(x, mu, sigma, nu){
  if (any(x <= 0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(mu <= 0)) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))  
  if (any(nu <= 0)) 
    stop(paste("nu must be positive", "\n", "")) 
  
  h <- dGIW(x, mu, sigma, nu, log=FALSE) / 
    pGIW(q=x, mu, sigma, nu, lower.tail=FALSE, log.p=FALSE)
  h 
}
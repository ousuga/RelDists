#' Additive Weibull distribution
#' 
#' @description 
#' Density, distribution function, quantile function, 
#' random generation and hazard function for the Additive Weibull distribution
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
#' Additive Weibull Distribution with parameters \code{mu}, 
#' \code{sigma}, \code{nu} and \code{tau} has density given by
#' 
#' \eqn{f(x) = (\mu\nu x^{\nu - 1} + \sigma\tau x^{\tau - 1}) \exp({-\mu x^{\nu} - \sigma x^{\tau} }),}
#' 
#' for x > 0. 
#' 
#' @return 
#' \code{dAdd} gives the density, \code{pAdd} gives the distribution 
#' function, \code{qAdd} gives the quantile function, \code{rAdd}
#' generates random deviates and \code{hAdd} gives the hazard function.
#'
#' @examples  
#' ## The probability density function
#' curve(dAdd(x, mu=1.5, sigma=0.5, nu=3, tau=0.8), from=0.0001, to=2,
#'       col="red", las=1, ylab="f(x)")
#' 
#' ## The cumulative distribution and the Reliability function
#' par(mfrow=c(1, 2))
#' curve(pAdd(x, mu=1.5, sigma=0.5, nu=3, tau=0.8),
#'       from=0.0001, to=2, col="red", las=1, ylab="F(x)")
#' curve(pAdd(x, mu=1.5, sigma=0.5, nu=3, tau=0.8, lower.tail=FALSE),
#'      from=0.0001, to=2, col="red", las=1, ylab="S(x)")
#' 
#' ## The quantile function
#' p <- seq(from=0, to=0.99999, length.out=100)
#' plot(x=qAdd(p, mu=1.5, sigma=0.2, nu=3, tau=0.8), y=p, xlab="Quantile",
#'      las=1, ylab="Probability")
#' curve(pAdd(x, mu=1.5, sigma=0.2, nu=3, tau=0.8), 
#'       from=0, add=TRUE, col="red")
#' 
#' ## The random function
#' hist(rAdd(n=10000, mu=1.5, sigma=0.2, nu=3, tau=0.8), freq=FALSE,
#'      xlab="x", las=1, main="")
#' curve(dAdd(x, mu=1.5, sigma=0.2, nu=3, tau=0.8),
#'       from=0.09, to=5, add=TRUE, col="red")
#' 
#' ## The Hazard function
#' curve(hAdd(x, mu=1.5, sigma=0.2, nu=3, tau=0.8), from=0.001, to=1,
#'       col="red", ylab="Hazard function", las=1)
#'
#' @export
dAdd <- function(x, mu, sigma,
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
  
  A <- mu*nu*x^(nu-1)
  B <- sigma*tau*x^(tau-1)
  loglik <- log(A+B) - (mu*x^(nu) + sigma*x^(tau))
  
  if (log == FALSE) 
    density <- exp(loglik)
  else 
    density <- loglik
  return(density)
}
#' @export
#' @rdname dAdd
pAdd <- function(q, mu, sigma, nu, tau, 
                 lower.tail=TRUE, log.p=FALSE){
  if (any(mu <= 0)) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", "")) 
  if (any(nu <= 0)) 
    stop(paste("nu must be positive", "\n", "")) 
  if (any(tau <= 0)) 
    stop(paste("tau must be postive", "\n", "")) 
  
  cdf <- 1 - exp(-mu*q^(nu) - sigma*q^(tau))
  
  if (lower.tail == TRUE) 
    cdf <- cdf
  else cdf <- 1 - cdf
  if (log.p == FALSE) 
    cdf <- cdf
  else cdf <- log(cdf)
  cdf
}
#' @export
#' @rdname dAdd
qAdd <- function(p, mu, sigma, nu, tau,
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
    1 - exp(-mu*x^(nu)-sigma*x^(tau))
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
#' @rdname dAdd
rAdd <- function(n, mu, sigma, nu, tau){
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
  r <- qAdd(p, mu, sigma, nu, tau)
  r
}
#' @export
#' @rdname dAdd
hAdd<-function(x, mu, sigma, nu, tau){
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
  
  h <- dAdd(x, mu, sigma, nu, tau, log=FALSE) / 
    pAdd(q=x, mu, sigma, nu, tau, lower.tail=FALSE, log.p=FALSE)
  h  
}


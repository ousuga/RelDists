#' The Four-Parameter Exponentiated Generalized Gamma distribution
#' 
#' @author Amylkar Urrea Montoya, \email{amylkar.urrea@@udea.edu.co}
#' 
#' @description 
#' Density, distribution function, quantile function, 
#' random generation and hazard function for the Four-Parameter Exponentiated Generalized Gamma distribution
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
#' \eqn{f(x) = \frac {\nu \sigma}{\mu \tau(\tau)} (\frac {x}{\mu})^{\sigma \tau -1} exp\{(-{\frac {x}{\mu})^{\sigma}})\} 
#' \{ \gamma~_1~ (\tau, (\frac {x}{\mu})^{\sigma}) \}^{\nu -1} ,}
#' 
#' for x > 0. 
#' 
#' @return 
#' \code{dFPEGG} gives the density, \code{pFPEGG} gives the distribution 
#' function, \code{qFPEGG} gives the quantile function, \code{rFPEGG}
#' generates random deviates and \code{hFPEGG} gives the hazard function.
#'
#' @examples  
#' ## The probability density function
#' curve(dFPEGG(x, mu=0.1, sigma=0.8, nu=10, tau=1.5), from=0.000001, to=1.5, ylim=c(0, 2.5),
#'       col="red", las=1, ylab="f(x)")
#' 
#' ## The cumulative distribution and the Reliability function
#' par(mfrow=c(1, 2))
#' curve(pFPEGG(x, mu=0.1, sigma=0.8, nu=10, tau=1.5),
#'       from=0.000001, to=1.5, col="red", las=1, ylab="F(x)")
#' curve(pFPEGG(x, mu=0.1, sigma=0.8, nu=10, tau=1.5, lower.tail=FALSE),
#'       from=0.000001, to=1.5, col="red", las=1, ylab="S(x)")
#' 
#' ## The quantile function
#' p <- seq(from=0, to=0.99999, length.out=100)
#' plot(x=qFPEGG(p, mu=0.1, sigma=0.8, nu=10, tau=1.5), y=p, xlab="Quantile",
#'      las=1, ylab="Probability")
#' curve(pFPEGG(x, mu=0.1, sigma=0.8, nu=10, tau=1.5), 
#'       from=0.00001, add=TRUE, col="red")
#' 
#' ## The random function
#' hist(rFPEGG(n=10000, mu=0.1, sigma=0.8, nu=10, tau=1.5), freq=FALSE,
#'      xlab="x", las=1, main="")
#' curve(dFPEGG(x, mu=0.1, sigma=0.8, nu=10, tau=1.5),
#'       from=0.0001, to=2, add=TRUE, col="red")
#' 
#' ## The Hazard function
#' curve(hFPEGG(x,  mu=0.1, sigma=0.8, nu=10, tau=1.5), from=0.0001, to=1.5,
#'       col="red", ylab="Hazard function", las=1)
#'
#' @references
#' \insertRef{almalki2014modifications}{RelDists}
#'
#' \insertRef{cordeiro2011}{RelDists}
#'
#' @importFrom Rdpack reprompt
#'
#' @export
dFPEGG <- function(x, mu, sigma,
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
#' @rdname dFPEGG
pFPEGG <- function(q, mu, sigma, nu, tau, 
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
#' @rdname dFPEGG
qFPEGG <- function(p, mu, sigma, nu, tau,
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
#' @rdname dFPEGG
rFPEGG <- function(n, mu, sigma, nu, tau){
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
  r <- qFPEGG(p, mu, sigma, nu, tau)
  r
}
#' @export
#' @rdname dFPEGG
hFPEGG <- function(x, mu, sigma, nu, tau){
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
  
  h <- dFPEGG(x, mu, sigma, nu, tau, log=FALSE) / 
    pFPEGG(q=x, mu, sigma, nu, tau, lower.tail=FALSE, log.p=FALSE)
  h  
}
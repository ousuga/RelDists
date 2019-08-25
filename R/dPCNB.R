#' The Power-Cauchy Negative-Binomial distribution
#' 
#' @author Santiago Uribe
#' 
#' @description 
#' Density, distribution function, quantile function, 
#' random generation and hazard function for the Power-Cauchy Negative-Binomial distribution
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
#' Power-Cauchy Negative-Binomial distribution with parameters \code{mu}, 
#' \code{sigma}, \code{nu} and \code{tau} has density given by
#' 
#' \eqn{f(x) = \frac{2\pi^{-1}\nu\tau^{\nu}(\mu/\sigma)(x/\sigma)^{\mu-1}\left[1+(x/\sigma)^{2\mu}\right]^{-1}
#' \left[2\pi^{-1}\tan^{-1}(x/\sigma)^{\mu}\right]^{\nu-1}}
#' {\left[1-2\pi^{-1}(1-\tau)\tan^{-1}(x/\sigma)^{\mu}\right]^{\nu+1}},}
#' 
#' for \eqn{x > 0}, \eqn{\mu> 0}, \eqn{\sigma> 0}, \eqn{\nu> 0}, \eqn{ 0 <\nu< 1}  
#' 
#' @return 
#' \code{dPCNB} gives the density, \code{pPCNB} gives the distribution 
#' function, \code{qPCNB} gives the quantile function, \code{rPCNB}
#' generates random deviates and \code{hPCNB} gives the hazard function.
#'
#' @examples  
#' ## The probability density function
#' curve(dPCNB(x, mu=3, sigma=1, nu=0.9, tau=0.8), from=0.001, 
#'       to=3, ylim=c(0, 1), col="red", las=1, ylab="f(x)")
#' 
#' ## The cumulative distribution and the Reliability function
#' par(mfrow=c(1, 2))
#' curve(pPCNB(x, mu=3, sigma=1, nu=0.9, tau=0.8),
#'       from=0.01, to=2.5, col="red", las=1, ylab="F(x)")
#' curve(pPCNB(x, mu=3, sigma=1, nu=0.9, tau=0.8, lower.tail=FALSE),
#'       from=0.01, to=2.5, col="red", las=1, ylab="S(x)")
#' 
#' ## The quantile function
#' p <- seq(from=0, to=0.99999, length.out=100)
#' plot(x=qPCNB(p, mu=3, sigma=1, nu=0.9, tau=0.8), y=p, xlab="Quantile",
#'      las=1, ylab="Probability")
#' curve(pPCNB(x, mu=3, sigma=1, nu=0.9, tau=0.8), 
#'       from=0.00001, add=TRUE, col="red")
#' 
#' ## The random function
#' hist(rPCNB(n=10000, mu=3, sigma=1, nu=0.9, tau=0.8), freq=FALSE,
#'      xlab="x", las=1, main="")
#' curve(dPCNB(x, mu=3, sigma=1, nu=0.9, tau=0.8),
#'       from=0.0001, to=10, add=TRUE, col="red")
#' 
#' ## The Hazard function
#' curve(hPCNB(x, mu=3, sigma=1, nu=0.9, tau=0.8), from=0.001, to=20,
#'       col="red", ylab="Hazard function", las=1)
#' 
#'
#' @references
#' \insertRef{muhammad2018}{RelDists}
#'
#' @importFrom Rdpack reprompt
#'
#' @export
dPCNB <- function(x, mu, sigma,
                  nu, tau, log=FALSE){
  if (any(x <= 0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(mu <= 0)) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", "")) 
  if (any(nu <= 0)) 
    stop(paste("nu must be positive", "\n", "")) 
  if (any(any(tau < 0) || any(tau > 1)))
    stop(paste("tau must be between 0 and 1", "\n", ""))
  
  A <- log(2/pi) + log(nu) + nu*log(tau) + log(mu/sigma) + (mu-1)*log(x/sigma)
  B <- -log(1+(x/sigma)^(2*mu)) + (nu-1)*log((2/pi)*atan((x/sigma)^(mu)))
  C <- -((nu+1)*log(1-((2/pi)*(1-tau)*atan((x/sigma)^(mu)))))
  loglik <- A + B + C
  
  if (log == FALSE) 
    density <- exp(loglik)
  else 
    density <- loglik
  return(density)
}
#' @export
#' @rdname dPCNB
pPCNB <- function(q, mu, sigma, nu, tau, 
                  lower.tail=TRUE, log.p=FALSE){
  if (any(mu <= 0)) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", "")) 
  if (any(nu <= 0)) 
    stop(paste("nu must be positive", "\n", "")) 
  if (any(any(tau < 0) || any(tau > 1)))
    stop(paste("tau must be between 0 and 1", "\n", ""))
  
  cdf <- (((2/pi)*tau*atan((q/sigma)^mu))/(1-((2/pi)*(1-tau)*atan((q/sigma)^mu))))^nu
  
  if (lower.tail == TRUE) 
    cdf <- cdf
  else cdf <- 1 - cdf
  if (log.p == FALSE) 
    cdf <- cdf
  else cdf <- log(cdf)
  cdf
}
#' @export
#' @rdname dPCNB
qPCNB <- function(p, mu, sigma, nu, tau,
                  lower.tail=TRUE, log.p=FALSE){
  if (any(mu <= 0)) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", "")) 
  if (any(nu <= 0)) 
    stop(paste("nu must be positive", "\n", "")) 
  if (any(any(tau < 0) || any(tau > 1)))
    stop(paste("tau must be between 0 and 1", "\n", ""))
  if (log.p == TRUE) 
    p <- exp(p)
  else p <- p
  if (lower.tail == TRUE) 
    p <- p
  else p <- 1 - p
  if (any(p < 0) | any(p > 1)) 
    stop(paste("p must be between 0 and 1", "\n", ""))
  
  fda <- function(x, mu, sigma, nu, tau){
    
    cdf <- (((2/pi)*tau*atan((x/sigma)^mu))/(1-((2/pi)*(1-tau)*atan((x/sigma)^mu))))^nu
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
#' @rdname dPCNB
rPCNB <- function(n, mu, sigma, nu, tau){
  if (any(mu <= 0)) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", "")) 
  if (any(nu <= 0)) 
    stop(paste("nu must be positive", "\n", "")) 
  if (any(any(tau < 0) || any(tau > 1)))
    stop(paste("tau must be between 0 and 1", "\n", ""))
  
  n <- ceiling(n)
  p <- runif(n)
  r <- qPCNB(p, mu, sigma, nu, tau)
  r
}
#' @export
#' @rdname dPCNB
hPCNB <- function(x, mu, sigma, nu, tau){
  if (any(x <= 0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(mu <= 0)) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", "")) 
  if (any(nu <= 0)) 
    stop(paste("nu must be positive", "\n", "")) 
  if (any(any(tau < 0) || any(tau > 1)))
    stop(paste("tau must be between 0 and 1", "\n", ""))
  
  h <- dPCNB(x, mu, sigma, nu, tau, log=FALSE) / 
    pPCNB(q=x, mu, sigma, nu, tau, lower.tail=FALSE, log.p=FALSE)
  h  
}

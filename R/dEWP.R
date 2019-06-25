#' The Exponentieted Weibull-Poisson distribution
#' 
#' @author Jaime Mosquera Gutiérrez, \email{jmosquerag@@unal.edu.co}
#' 
#' @description 
#' Density, distribution function, quantile function, 
#' random generation and hazard function for the odd Weibull distribution with
#' parameters \code{mu}, \code{sigma}, \code{nu} and \code{tau}.
#' 
#' @param x,q	vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. 
#' @param mu scale parameter.
#' @param sigma,nu,tau shape parameters.
#' @param log,log.p	logical; if TRUE, probabilities p are given as log(p).	
#' @param lower.tail logical; if TRUE (default), probabilities are P[X <= x], otherwise, P[X > x].
#' 
#' @details 
#' The exponentieted Weibull-Poisson distribution with parameters \code{mu}, 
#' \code{sigma}, \code{nu} and \code{tau} has density given by
#' 
#' \eqn{f(x)=\sigma \nu \tau \mu^\nu x^(\nu - 1) \exp(-(\mu x)^\nu) (1 - \exp(-(\mu x)^\nu))^(\sigma - 1)
#'            \exp(\tau(1 - \exp(-(\mu x)^\nu))^\sigma),}
#'            
#' for \eqn{x > 0}, \eqn{\mu > 0}, \eqn{\sigma > 0}, \eqn{\nu >0} and \eqn{\nu >0}.
#'
#' @return 
#' \code{dEWP} gives the density, \code{pEWP} gives the distribution 
#' function, \code{qEWP} gives the quantile function, \code{rEWP}
#' generates random deviates and \code{hEWP} gives the hazard function.
#' 
#' @examples  
#' ## The probability density function
#' par(mfrow = c(1, 2))
#' curve(dEWP(x, mu=1, sigma=0.5, nu=1, tau=5), from=0, to=4, ylim=c(0, 0.6), 
#'       col="red", las=1, ylab="The probability density function")
#' curve(dEWP(x, mu=1, sigma=0.5, nu=3, tau=0.5), from=0, to=2, ylim=c(0, 0.25), 
#'       col="red", las=1, ylab="The probability density function")
#'       
#' ## The cumulative distribution and the Reliability function
#' par(mfrow = c(1, 2))
#' curve(pEWP(x, mu=1, sigma=0.5, nu=1, tau=5), from=0, to=4, ylim=c(0, 1), 
#'       col="red", las = 1, ylab = "The cumulative distribution function")
#' curve(pEWP(x, mu=1, sigma=0.5, nu=1, tau=5, lower.tail=FALSE), from=0, to=10,  
#'       ylim=c(0, 1), col="red", las=1, ylab="The Reliability function")
#' 
#' @references
#' \insertRef{Mahmoudi2013}{RelDists}
#'
#' @importFrom Rdpack reprompt
#' 
#' @export
dEWP <- function(x, mu, sigma, nu, tau, log=FALSE){
  if (any(x < 0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(mu <= 0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma*nu <= 0)) 
    stop(paste("Product sigma*nu must be positive", "\n", ""))
  
  loglik <- log(sigma) + log(nu) + nu*log(mu) + (nu - 1)*log(x) - 
    log(exp(tau) - 1) - (mu*x)^nu + (sigma - 1)*log(1 - exp(-(mu*x)^nu)) +
    tau*(1 - exp(-(mu*x)^nu))^sigma
  
  if (log == FALSE)
    density<- exp(loglik)
  else 
    density <- loglik
  return(density)
}
#' @export
#' @rdname dEWP
pEWP <- function(q, mu, sigma, nu, tau, lower.tail=TRUE, log.p=FALSE){
  if (any(q < 0)) 
    stop(paste("q must be positive", "\n", ""))
  if (any(mu <= 0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma*nu <= 0)) 
    stop(paste("Product sigma*nu must be positive", "\n", ""))
  
  cdf <- exp(tau*(1 - exp((1 - exp(-(mu*q)^nu))^sigma)))/( exp(tau) - 1)
  if (lower.tail == TRUE)
    cdf <- cdf
  else cdf <- 1 - cdf
  if (log.p == FALSE) 
    cdf <- cdf
  else cdf <- log(cdf)
  cdf
}
#' @export
#' @rdname dEWP
qEWP <- function(p, mu, sigma, nu, tau, lower.tail=TRUE, log.p=FALSE){
  if (any(mu <= 0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma*nu <= 0)) 
    stop(paste("Product sigma*nu must be positive", "\n", ""))
  
  if (log.p == TRUE) 
    p <- exp(p)
  else p <- p
  if (lower.tail) 
    p <- p
  else p <- 1 - p
  if (any(p < 0) | any(p > 1)) 
    stop(paste("p must be between 0 and 1", "\n", ""))
  
  q <- (1/mu)*(-log(1 - ((1/tau)*log(p*(exp(tau) - 1) + 1))^(1/sigma)))^(1/nu)
  q
}
#' @importFrom stats runif
#' @export
#' @rdname dEWP
rEWP <- function(n, mu, sigma, nu, tau){
if(any(n <= 0))
  stop(paste("n must be a positive integer","\n",""))
if (any(mu <= 0 )) 
  stop(paste("mu must be positive", "\n", ""))
if (any(sigma*nu <= 0)) 
  stop(paste("Product sigma*nu must be positive", "\n", ""))

n <- ceiling(n)
p <- runif(n)
r <- qEWP(p, mu, sigma, nu, tau)
r
}
#' @export
#' @rdname dEWP
hEWP <- function(x, mu, sigma, nu, tau){
  if (any(x < 0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(mu <= 0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma*nu <= 0)) 
    stop(paste("Product sigma*nu must be positive", "\n", ""))
  
  h <- dEWP(x, mu, sigma, nu, tau, log=FALSE)/
    pEWP(q=x, mu, sigma, nu, tau, lower.tail=FALSE, log.p=FALSE)
  h
}
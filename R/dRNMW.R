#' The Reduced New Modified Weibull distribution
#' 
#' @author Amylkar Urrea Montoya, \email{amylkar.urrea@@udea.edu.co}
#' 
#' @description 
#' Density, distribution function, quantile function, 
#' random generation and hazard function for the Reduced New Modified Weibull distribution
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
#' The Reduced New Modified Weibull distribution \code{mu}, 
#' \code{sigma} and \code{nu} has density given by
#' 
#' \eqn{f(x) = \frac{1}{2 \sqrt{x}} \{\mu + \sigma(1+2 \nu x)\exp(\nu x) \}exp(-\mu \sqrt{x}-\sigma \sqrt{x} \exp(\nu x)),}
#' 
#' for x > 0. 
#' 
#' @return 
#' \code{dRNMW} gives the density, \code{pRNMW} gives the distribution 
#' function, \code{qRNMW} gives the quantile function, \code{rRNMW}
#' generates random deviates and \code{hRNMW} gives the hazard function.
#'
#' @examples  
#' ## The probability density function
#' curve(dRNMW(x, mu=0.5, sigma=0.015, nu=2.75), from=0.001, to=5,
#'       col="red", ylab="f(x)", las=1, ylim=c(0, 1))
#' 
#' ## The cumulative distribution and the Reliability function
#' par(mfrow=c(1, 2))
#' curve(pRNMW(x, mu=0.5, sigma=0.015, nu=2.751),
#'       from=0.01, to=2.5, col="red", las=1, ylab="F(x)")
#' curve(pRNMW(x, mu=0.5, sigma=0.015, nu=2.75, lower.tail=FALSE),
#'       from=0.01, to=2.5, col="red", las=1, ylab="S(x)")
#' 
#' ## The quantile function
#' p <- seq(from=0, to=0.99999, length.out=100)
#' plot(x=qRNMW(p, mu=0.5, sigma=0.015, nu=2.75), y=p, xlab="Quantile",
#'      las=1, ylab="Probability")
#' curve(pRNMW(x, mu=0.5, sigma=0.015, nu=2.75),
#'       from=0, add=TRUE, col="red")
#' 
#' ## The random function
#' hist(rRNMW(n=10000, mu=0.5, sigma=0.015, nu=2.75), freq=FALSE,
#'      xlab="x", ylim=c(0, 2), las=1, main="")
#' curve(dRNMW(x, mu=0.5, sigma=0.015, nu=2.75),
#'       from=0.001, to=4, add=TRUE, col="red")
#' 
#' ## The Hazard function
#' curve(hRNMW(x, mu=0.5, sigma=0.015, nu=2.75), from=0.001, to=3,
#'       col="red", ylab="Hazard function", las=1)
#'
#' @references
#' \insertRef{almalki2014modifications}{RelDists}
#' 
#' \insertRef{almalki2018}{RelDists}
#'
#' @importFrom Rdpack reprompt
#'
#' @export
dRNMW <- function(x, mu, sigma, nu, log=FALSE){
  if (any(x <= 0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(mu <= 0)) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", "")) 
  if (any(nu <= 0)) 
    stop(paste("nu must be positive", "\n", "")) 
  
  A <- log(1) - log(2 * sqrt(x)) - mu * sqrt(x)
  B <- -sigma * sqrt(x) * exp(nu * x)
  C <- log(mu + sigma * (1 + 2 * nu * x) * exp(nu * x))
  loglik <- A + B + C
  
  if (log == FALSE) 
    density <- exp(loglik)
  else 
    density <- loglik
  return(density)
}
#' @export
#' @rdname dRNMW
pRNMW <- function(q, mu, sigma, nu, 
                  lower.tail=TRUE, log.p=FALSE){
  if (any(mu <= 0)) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", "")) 
  if (any(nu <= 0)) 
    stop(paste("nu must be positive", "\n", "")) 
  
  A <- exp(-mu * sqrt(q) - sigma * sqrt(q) * exp(nu * q))
  cdf <- 1 - A
  
  if (lower.tail == TRUE) 
    cdf <- cdf
  else cdf <- 1 - cdf
  if (log.p == FALSE) 
    cdf <- cdf
  else cdf <- log(cdf)
  cdf
}
#' @export
#' @rdname dRNMW
qRNMW <- function(p, mu, sigma, nu,
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
    
    1 - exp(-mu * sqrt(x) - sigma * sqrt(x) * exp(nu * x))
    
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
#' @rdname dRNMW
rRNMW <- function(n, mu, sigma, nu){
  if (any(mu <= 0)) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu <= 0)) 
    stop(paste("nu must be positive", "\n", "")) 
  
  n <- ceiling(n)
  p <- runif(n)
  r <- qRNMW(p, mu, sigma, nu)
  r
}
#' @export
#' @rdname dRNMW
hRNMW<-function(x, mu, sigma, nu){
  if (any(x <= 0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(mu <= 0)) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))  
  if (any(nu <= 0)) 
    stop(paste("nu must be positive", "\n", "")) 
  
  h <- dRNMW(x, mu, sigma, nu, log=FALSE) / 
    pRNMW(q=x, mu, sigma, nu, lower.tail=FALSE, log.p=FALSE)
  h 
}

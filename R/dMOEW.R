#' The Marshall-Olkin Extended Weibull distribution
#' 
#' @author Amylkar Urrea Montoya, \email{amylkar.urrea@@udea.edu.co}
#' 
#' @description 
#' Density, distribution function, quantile function, 
#' random generation and hazard function for the Marshall-Olkin Extended Weibull distribution
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
#' The Marshall-Olkin Extended Weibull distribution \code{mu}, 
#' \code{sigma} and \code{nu} has density given by
#' 
#' \eqn{f(x) = \frac{\mu \sigma \nu (\nu x)^{\sigma - 1} exp\{{-(\nu x )^{\sigma}}\}}{\{1-(1-\mu) exp\{{-(\nu x )^{\sigma}}\} \}^{2}},}
#' 
#' for x > 0. 
#' 
#' @return 
#' \code{dMOEW} gives the density, \code{pMOEW} gives the distribution 
#' function, \code{qMOEW} gives the quantile function, \code{rMOEW}
#' generates random deviates and \code{hMOEW} gives the hazard function.
#'
#' @examples  
#' ## The probability density function
#' curve(dMOEW(x, mu=0.5, sigma=0.7, nu=1), from=0.001, to=10,
#'       col="red", ylab="f(x)", las=1)
#' 
#' ## The cumulative distribution and the Reliability function
#' par(mfrow=c(1, 2))
#' curve(pMOEW(x, mu=0.5, sigma=0.7, nu=1),
#'       from=0.0001, to=2, col="red", las=1, ylab="F(x)")
#' curve(pMOEW(x, mu=0.5, sigma=0.7, nu=1, lower.tail=FALSE),
#'       from=0.0001, to=2, col="red", las=1, ylab="S(x)")
#' 
#' ## The quantile function
#' p <- seq(from=0, to=0.99999, length.out=100)
#' plot(x=qMOEW(p, mu=0.5, sigma=0.7, nu=1), y=p, xlab="Quantile",
#'      las=1, ylab="Probability")
#' curve(pMOEW(x, mu=0.5, sigma=0.7, nu=1),
#'       from=0, add=TRUE, col="red")
#' 
#' ## The random function
#' hist(rMOEW(n=10000, mu=0.5, sigma=0.7, nu=1), freq=FALSE,
#'      xlab="x", ylim=c(0, 2), las=1, main="")
#' curve(dMOEW(x, mu=0.5, sigma=0.7, nu=1),
#'       from=0.001, to=4, add=TRUE, col="red")
#' 
#' ## The Hazard function
#' curve(hMOEW(x, mu=0.5, sigma=0.7, nu=1), from=0.001, to=3,
#'       col="red", ylab="Hazard function", las=1)
#'
#' @references
#' \insertRef{almalki2014modifications}{RelDists}
#' 
#' \insertRef{ghitany2005}{RelDists}
#'
#' @importFrom Rdpack reprompt
#'
#' @export
dMOEW <- function(x, mu, sigma, nu, log=FALSE){
  if (any(x <= 0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(mu <= 0)) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", "")) 
  if (any(nu <= 0)) 
    stop(paste("nu must be positive", "\n", "")) 
  
  A <- log(mu) + log(sigma) + log(nu) + (sigma - 1) * log(nu * x) - (nu * x)^sigma  
  B <- 2 * log(1 - (1 - mu) * exp(-(nu * x)^sigma))
  loglik <- A - B
  
  if (log == FALSE) 
    density <- exp(loglik)
  else 
    density <- loglik
  return(density)
}
#' @export
#' @rdname dMOEW
pMOEW <- function(q, mu, sigma, nu, 
                  lower.tail=TRUE, log.p=FALSE){
  if (any(mu <= 0)) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", "")) 
  if (any(nu <= 0)) 
    stop(paste("nu must be positive", "\n", "")) 
  
  A <- mu * exp(-(nu * q)^sigma)
  B <- 1 - (1 - mu) * exp(-(nu * q)^sigma)
  cdf <- 1 - (A / B)
  
  if (lower.tail == TRUE) 
    cdf <- cdf
  else cdf <- 1 - cdf
  if (log.p == FALSE) 
    cdf <- cdf
  else cdf <- log(cdf)
  cdf
}
#' @export
#' @rdname dMOEW
qMOEW <- function(p, mu, sigma, nu,
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
    
    1 - (mu * exp(-(nu * x)^sigma)) / (1 - (1 - mu) * exp(-(nu * x)^sigma)) 
    
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
#' @rdname dMOEW
rMOEW <- function(n, mu, sigma, nu){
  if (any(mu <= 0)) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu <= 0)) 
    stop(paste("nu must be positive", "\n", "")) 
  
  n <- ceiling(n)
  p <- runif(n)
  r <- qMOEW(p, mu, sigma, nu)
  r
}
#' @export
#' @rdname dMOEW
hMOEW<-function(x, mu, sigma, nu){
  if (any(x <= 0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(mu <= 0)) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))  
  if (any(nu <= 0)) 
    stop(paste("nu must be positive", "\n", "")) 
  
  h <- dMOEW(x, mu, sigma, nu, log=FALSE) / 
    pMOEW(q=x, mu, sigma, nu, lower.tail=FALSE, log.p=FALSE)
  h 
}

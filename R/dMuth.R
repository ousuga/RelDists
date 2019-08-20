#' The Muth distribution
#' 
#' @author William Aguilar
#' 
#' @description 
#' Density, distribution function, quantile function, 
#' random generation and hazard function for the Power Lindley distribution 
#' with parameter \code{mu}.
#' 
#' @param x,q	vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. 
#' @param mu parameter.
#' @param log,log.p	logical; if TRUE, probabilities p are given as log(p).	
#' @param lower.tail logical; if TRUE (default), probabilities are P[X <= x], otherwise, P[X > x].
#' 
#' @details 
#' The Muth distribution with parameter \code{mu} 
#' has density given by
#' 
#' \eqn{f(x)= ({exp(\mu x) - \mu})exp[{\mu x - \frac{1}{\mu} ({exp(\mu x) - 1})}]}
#' 
#' for \eqn{x > 0}, \eqn{0< \mu <=1}.
#' 
#' @return 
#' \code{dMuth} gives the density, \code{pMuth} gives the distribution 
#' function, \code{qMuth} gives the quantile function, \code{rMuth}
#' generates random deviates and \code{hMuth} gives the hazard function.
#'
#' @examples  
#' ## The probability density function
#' curve(dMuth(x, mu=1), from=0.1, to=5,
#'       ylim=c(0, 1), col="red", las=1, ylab="f(x)")
#' 
#' ## The cumulative distribution and the Reliability function
#' par(mfrow=c(1, 2))
#' curve(pMuth(x, mu=1),
#'       from=0.1, to=2,  col="red", las=1, ylab="F(x)")
#' curve(pMuth(x,  mu=1, lower.tail=FALSE),
#'       from=0.1, to=2, col="red", las=1, ylab="S(x)")
#' 
#' ## The quantile function
#' p <- seq(from=0, to=0.99999, length.out=100)
#' plot(x=qMuth(p, mu=1), y=p, xlab="Quantile",
#'      las=1, ylab="Probability")
#' curve(pMuth(x, mu=1), from=0.01, add=TRUE, col="red")
#' 
#' ## The random function
#' hist(rMuth(n=10000, mu=1), freq=FALSE,
#'      xlab="x", las=1, main="")
#' curve(dMuth(x, mu=1), from=0.1, to=8, add=TRUE, col="red")
#' 
#' ## The Hazard function
#' par(mfrow=c(1, 1))
#' curve(hMuth(x, mu=1), from=0.1, to=3.5,
#'       col="red", ylab="Hazard function", las=1)
#'
#' @references
#' \insertRef{abdullah2018}{RelDists}
#'
#' @importFrom Rdpack reprompt
#'
#' @export
dMuth <- function(x, mu, log=FALSE){
  if (any(x <= 0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(mu <= 0) | any(mu > 1)) 
    stop(paste("mu must be between into (0,1]", "\n", ""))
  
  loglik <- log(exp(mu * x) -mu) + (mu * x - (1/mu) * (exp(mu * x)-1))
  
  if (log == FALSE)
    density <- exp(loglik)
  else
    density <- loglik
  return(density)
}
#' @export
#' @rdname dMuth
pMuth <- function(q, mu, lower.tail=TRUE, log.p=FALSE){
  if (any(q <= 0)) 
    stop(paste("q must be positive", "\n", ""))
  if (any(mu <= 0) | any(mu > 1)) 
    stop(paste("mu must be between into (0,1]", "\n", ""))
  
  cdf <- 1 - exp(mu * q - (1/mu) * (exp(mu * q) - 1))
  
  if (lower.tail == TRUE) 
    cdf <- cdf
  else cdf <- 1 - cdf
  if (log.p == FALSE) 
    cdf <- cdf
  else cdf <- log(cdf)
  cdf
}
#' @export
#' @rdname dMuth
qMuth <- function(p, mu, lower.tail=TRUE, log.p=FALSE){
  if (any(mu <= 0) | any(mu > 1)) 
    stop(paste("mu must be between into (0,1]", "\n", ""))
  
  if (log.p == TRUE) 
    p <- exp(p)
  else p <- p
  if (lower.tail == TRUE) 
    p <- p
  else p <- 1 - p
  if (any(p < 0) | any(p > 1)) 
    stop(paste("p must be between 0 and 1", "\n", ""))
  
  fda <- function(x, mu){
    
    1 - exp(mu * x - (1/mu) * (exp(mu * x) - 1))
    
  }
  fda1 <- function(x, mu, p) {
    fda(x, mu) - p
  }
  r_de_la_funcion <- function(mu, p) {
    uniroot(fda1, interval=c(0, 1e+06), mu, p)$root
  }
  r_de_la_funcion <- Vectorize(r_de_la_funcion)
  q <- r_de_la_funcion(mu, p)
  q
}
#' @importFrom stats runif
#' @export
#' @rdname dMuth
rMuth <- function(n, mu, log.p=FALSE){
  if(any(n <= 0))
    stop(paste("n must be positive","\n",""))
  if (any(mu <= 0) | any(mu > 1)) 
    stop(paste("mu must be between into (0,1]", "\n", ""))
  
  n <- ceiling(n)
  p <- runif(n)
  r <- qMuth(p, mu)
  r
}
#' @export
#' @rdname dMuth
hMuth<-function(x, mu){
  if (any(x <= 0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(mu <= 0) | any(mu > 1)) 
    stop(paste("mu must be between into (0,1]", "\n", ""))
  
  h <- dMuth(x, mu, log=FALSE) / 
    pMuth(q=x, mu, lower.tail=FALSE, log.p=FALSE)
  h  
}
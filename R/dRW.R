#' The Reflected Weibull distribution
#' 
#' @author Amylkar Urrea Montoya, \email{amylkar.urrea@@udea.edu.co}
#' 
#' @description 
#' Density, distribution function, quantile function, 
#' random generation and hazard function for the Reflected Weibull Distribution 
#' with parameters \code{mu} and \code{sigma}.
#' 
#' @param x,q	vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. 
#' @param mu parameter.
#' @param sigma parameter.
#' @param log,log.p	logical; if TRUE, probabilities p are given as log(p).	
#' @param lower.tail logical; if TRUE (default), probabilities are P[X <= x], otherwise, P[X > x].
#' 
#' @details 
#' The Reflected Weibull Distribution with parameters \code{mu} 
#' and \code{sigma} has density given by
#' 
#' \eqn{f(y) = \mu\sigma (-y) ^{\sigma - 1} e ^ {-\mu(-y)^\sigma},}
#' 
#' for y < 0.
#' 
#' @return 
#' \code{dRW} gives the density, \code{pRW} gives the distribution 
#' function, \code{qRW} gives the quantile function, \code{rRW}
#' generates random deviates and \code{hRW} gives the hazard function.
#'
#' @examples  
#' ## The probability density function
#' curve(dRW(x, mu=1, sigma=1), from=-5, to=-0.01,
#'       col="red", las=1, ylab="f(x)")
#' 
#' ## The cumulative distribution and the Reliability function
#' par(mfrow=c(1, 2))
#' curve(pRW(x, mu=1, sigma=1),
#'       from=-5, to=-0.01, col="red", las=1, ylab="F(x)")
#' curve(pRW(x, mu=1, sigma=1, lower.tail=FALSE),
#'       from=-5, to=-0.01, col="red", las=1, ylab="S(x)")
#' 
#' ## The quantile function
#' p <- seq(from=0, to=0.99999, length.out=100)
#' plot(x=qRW(p, mu=1, sigma=1), y=p, xlab="Quantile",
#'      las=1, ylab="Probability")
#' curve(pRW(x, mu=1, sigma=1), from=-5, add=TRUE, col="red")
#' 
#' ## The random function
#' hist(rRW(n=10000, mu=1, sigma=1), freq=FALSE,
#'      xlab="x", las=1, main="")
#' curve(dRW(x, mu=1, sigma=1), from=-5, to=-0.01, add=TRUE, col="red")
#' 
#' ## The Hazard function
#' par(mfrow=c(1,1))
#' curve(hRW(x, mu=1, sigma=1), from=-0.3, to=-0.01,
#'       col="red", ylab="Hazard function", las=1)
#'
#' @references
#' \insertRef{almalki2014modifications}{RelDists}
#' 
#' \insertRef{Clifford1973}{RelDists}
#'
#' @importFrom Rdpack reprompt
#'
#' @export
dRW <- function(x, mu, sigma, log=FALSE){
  if (any(x >= 0)) 
    stop(paste("x must be negative", "\n", ""))
  if (any(mu <= 0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  
  loglik<- log(mu) + log(sigma) + (sigma-1)*log(-x) -
    mu*((-x)^sigma)
  
  if (log == FALSE) 
    density <- exp(loglik)
  else 
    density <- loglik
  return(density)
}
#' @export
#' @rdname dRW
pRW <- function(q, mu, sigma,
                lower.tail=TRUE, log.p=FALSE){
  if (any(q >= 0)) 
    stop(paste("q must be negative", "\n", ""))
  if (any(mu <= 0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  
  cdf <- exp(-mu*(-q)^sigma)
  if (lower.tail == TRUE) 
    cdf <- cdf
  else cdf <- 1 - cdf
  if (log.p == FALSE) 
    cdf <- cdf
  else cdf <- log(cdf)
  cdf
}
#' @export
#' @rdname dRW
qRW <- function(p, mu, sigma,
                lower.tail=TRUE, log.p=FALSE){
  if (any(mu <= 0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  
  if (log.p == TRUE) 
    p <- exp(p)
  else p <- p
  if (lower.tail == TRUE) 
    p <- p
  else p <- 1 - p
  if (any(p < 0) | any(p > 1)) 
    stop(paste("p must be between 0 and 1", "\n", ""))
  
  q <- -{((-1/mu)*log(p))^(1/sigma)}
  q
}
#' @importFrom stats runif
#' @export
#' @rdname dRW
rRW <- function(n, mu, sigma){
  if(any(n <= 0))
    stop(paste("n must be positive","\n",""))
  if (any(mu <= 0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  
  n <- ceiling(n)
  p <- runif(n)
  r <- qRW(p, mu, sigma)
  r
}
#' @export
#' @rdname dRW
hRW<-function(x, mu, sigma){
  if (any(x >= 0)) 
    stop(paste("x must be negative", "\n", ""))
  if (any(mu <= 0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  
  h <- dRW(x, mu, sigma, log=FALSE) /
    pRW(q=x, mu, sigma, lower.tail=FALSE, log.p=FALSE)
  h
}


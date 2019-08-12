#' The Log-Weibull distribution
#' 
#' @author Amylkar Urrea Montoya, \email{amylkar.urrea@@udea.edu.co}
#' 
#' @description 
#' Density, distribution function, quantile function, 
#' random generation and hazard function for the Log-Weibull distribution 
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
#' The Log-Weibull Distribution with parameters \code{mu} 
#' and \code{sigma} has density given by
#' 
#' \eqn{f(y)=(1/\sigma) e^{((y - \mu)/\sigma)} exp\{-e^{((y - \mu)/\sigma)}\},}
#' 
#' for - \code{infty} < y < \code{infty}. 
#' 
#' @return 
#' \code{dLW} gives the density, \code{pLW} gives the distribution 
#' function, \code{qLW} gives the quantile function, \code{rLW}
#' generates random deviates and \code{hLW} gives the hazard function.
#'
#' @examples  
#' ## The probability density function
#' curve(dLW(x, mu=0, sigma=1.5), from=-8, to=5,
#'       col="red", las=1, ylab="f(x)")
#' 
#' ## The cumulative distribution and the Reliability function
#' par(mfrow=c(1, 2))
#' curve(pLW(x, mu=0, sigma=1.5),
#'       from=-8, to=5,  col="red", las=1, ylab="F(x)")
#' curve(pLW(x, mu=0, sigma=1.5, lower.tail=FALSE),
#'       from=-8, to=5, col="red", las=1, ylab="S(x)")
#' 
#' ## The quantile function
#' p <- seq(from=0, to=0.99999, length.out=100)
#' plot(x=qLW(p, mu=0, sigma=1.5), y=p, xlab="Quantile",
#'      las=1, ylab="Probability")
#' curve(pLW(x, mu=0, sigma=1.5), from=-8, to=5, add=TRUE, col="red")
#' 
#' ## The random function
#' hist(rLW(n=10000, mu=0, sigma=1.5), freq=FALSE,
#'      xlab="x", las=1, main="")
#' curve(dLW(x, mu=0, sigma=1.5), from=-15, to=6, add=TRUE, col="red")
#' 
#' ## The Hazard function
#' par(mfrow=c(1,1))
#' curve(hLW(x, mu=0, sigma=1.5), from=-8, to=7,
#'       col="red", ylab="Hazard function", las=1)
#'
#' @references
#' \insertRef{almalki2014modifications}{RelDists}
#' 
#' \insertRef{Gumbel1958}{RelDists}
#'
#' @importFrom Rdpack reprompt
#'
#' @export
dLW <- function(x, mu, sigma, log=FALSE){
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  
  loglik <- -log(sigma) + (x-mu)/sigma - exp((x-mu)/sigma)
  
  if (log == FALSE) 
    density <- exp(loglik)
  else density <- loglik
  return(density)
}
#' @export
#' @rdname dLW
pLW <- function(q, mu, sigma, 
                lower.tail=TRUE, log.p=FALSE){
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  
  cdf <- 1 - exp(-exp((q-mu)/sigma))
  
  if (lower.tail == TRUE) 
    cdf <- cdf
  else cdf <- 1 - cdf
  if (log.p == FALSE) 
    cdf <- cdf
  else cdf <- log(cdf)
  cdf
}
#' @export
#' @rdname dLW
qLW <- function(p, mu, sigma, lower.tail=TRUE, log.p=FALSE){
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
  
  q <- sigma*(log(-log(1-p))) + mu  
  q
}
#' @importFrom stats runif
#' @export
#' @rdname dLW
rLW <- function(n, mu, sigma){
  if(any(n <= 0))
    stop(paste("n must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  
  n <- ceiling(n)
  p <- runif(n)
  r <- qLW(p, mu, sigma)
  r
}
#' @export
#' @rdname dLW
hLW<-function(x, mu, sigma){
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  
  h <- dLW(x, mu, sigma, log = FALSE) / 
    pLW(q=x, mu, sigma, lower.tail=FALSE, log.p=FALSE)
  h
}
#' The Power Lindley distribution
#' 
#' @description 
#' Density, distribution function, quantile function, 
#' random generation and hazard function for the Power Lindley distribution 
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
#' The Power Lindley Distribution with parameters \code{mu} 
#' and \code{sigma} has density given by
#' 
#' \eqn{f(x) = \frac{\mu \sigma^2}{\sigma + 1} (1 + x^\mu) x ^ {\mu - 1} \exp({-\sigma x ^\mu}),}
#' 
#' for x > 0.
#' 
#' @return 
#' \code{dPL} gives the density, \code{pPL} gives the distribution 
#' function, \code{qPL} gives the quantile function, \code{rPL}
#' generates random deviates and \code{hPL} gives the hazard function.
#'
#' @examples  
#' ## The probability density function
#' curve(dPL(x, mu=1.5, sigma=0.2), from=0.1, to=10,
#'       col="red", las=1, ylab="f(x)")
#' 
#' ## The cumulative distribution and the Reliability function
#' par(mfrow=c(1, 2))
#' curve(pPL(x, mu=1.5, sigma=0.2),
#'       from=0.1, to=10, col="red", las=1, ylab="F(x)")
#' curve(pPL(x, mu=1.5, sigma=0.2, lower.tail=FALSE),
#'       from=0.1, to=10, col="red", las=1, ylab="S(x)")
#' 
#' ## The quantile function
#' p <- seq(from=0, to=0.99999, length.out=100)
#' plot(x=qPL(p, mu=1.5, sigma=0.2), y=p, xlab="Quantile",
#'      las=1, ylab="Probability")
#' curve(pPL(x, mu=1.5, sigma=0.2), from=0.1, add=TRUE, col="red")
#' 
#' ## The random function
#' hist(rPL(n=10000, mu=1.5, sigma=0.2), freq=FALSE,
#'      xlab="x", las=1, main="")
#' curve(dPL(x, mu=1.5, sigma=0.2), from=0.1, to=15, add=TRUE, col="red")
#' 
#' ## The Hazard function
#' par(mfrow=c(1,1))
#' curve(hPL(x, mu=1.5, sigma=0.2), from=0.1, to=15,
#'       col="red", ylab="Hazard function", las=1)
#'
#' @export
dPL <- function(x, mu, sigma, log=FALSE){
  if (any(x <= 0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(mu <= 0)) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", "")) 
  
  loglik <- log(mu) + 2*log(sigma) - log(sigma+1) +
    log(1+(x^mu)) + (mu-1)*log(x) - sigma*(x^mu)
  
  if (log == FALSE) 
    density <- exp(loglik)
  else density <- loglik
  return(density)
}
#' @export
#' @rdname dPL
pPL <- function(q, mu, sigma, 
                lower.tail=TRUE, log.p=FALSE){
  if (any(q <= 0)) 
    stop(paste("q must be positive", "\n", ""))
  if (any(mu <= 0)) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))  
  
  cdf <- 1 - (1+((sigma/(sigma+1))*q^mu))*exp(-sigma*(q^mu))
  
  if (lower.tail == TRUE) 
    cdf <- cdf
  else cdf <- 1 - cdf
  if (log.p == FALSE) 
    cdf <- cdf
  else cdf <- log(cdf)
  cdf
}
#' @export
#' @rdname dPL
qPL <- function(p, mu, sigma, 
                lower.tail=TRUE, log.p=FALSE){
  if (any(mu <= 0)) 
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
  
  fda <- function(x, mu, sigma){
    1 - (1+((sigma/(sigma+1))*x^mu))*exp(-sigma*(x^mu))
  }
  fda1 <- function(x, mu, sigma, p) {
    fda(x, mu, sigma) - p
  }
  r_de_la_funcion <- function(mu, sigma, p) {
    uniroot(fda1, interval=c(0, 1e+06), mu, sigma, p)$root
  }
  r_de_la_funcion <- Vectorize(r_de_la_funcion)
  q <- r_de_la_funcion(mu, sigma, p)
  q
}
#' @importFrom stats runif
#' @export
#' @rdname dPL
rPL <- function(n, mu, sigma){
  if(any(n <= 0))
    stop(paste("n must be positive","\n",""))
  if (any(mu <= 0)) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  
  n <- ceiling(n)
  p <- runif(n)
  r <- qPL(p, mu, sigma)
  r
}
#' @export
#' @rdname dPL
hPL<-function(x, mu, sigma){
  if (any(x <= 0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(mu <= 0)) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  
  h <- dPL(x, mu, sigma, log=FALSE) / 
    pPL(q=x, mu, sigma, lower.tail=FALSE, log.p=FALSE)
  h  
}
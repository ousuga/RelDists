#' The Alpha Power Inverted Exponential distribution
#' 
#' @author Santiago Restrepo Alvarez
#' 
#' @description 
#' Density, distribution function, quantile function, 
#' random generation and hazard function for the Alpha Power Inverted Exponential distribution
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
#' The Alpha Power Inverted Exponential distribution \code{mu} and 
#' \code{sigma}  has density given by
#' 
#' \eqn{f(x) = \frac{log(\mu)}{\mu - 1} \frac{\sigma}{x^{2}} \exp(-\frac{\sigma}{x}) 
#' \mu^{\exp(-\frac{\sigma}{x})},}
#' 
#' for \eqn{x > 0}, \eqn{\mu> 0}, \eqn{\sigma> 0}, \eqn{\mu \neq 1}.
#' 
#' \eqn{f(x) = \frac{\sigma}{x^{2}} \exp(-\frac{\sigma}{x}),}
#' 
#' for \eqn{x > 0}, \eqn{\mu> 0}, \eqn{\sigma> 0}, \eqn{\mu = 1}.
#' 
#' @return 
#' \code{dAPIE} gives the density, \code{pAPIE} gives the distribution 
#' function, \code{qAPIE} gives the quantile function, \code{rAPIE}
#' generates random deviates and \code{hAPIE} gives the hazard function.
#'
#' @examples  
#' ## The probability density function
#' curve(dAPIE(x, mu=3, sigma=0.5), from=0.001, to=10,
#'       ylim=c(0, 0.8), col="red", ylab="f(x)", las=1)
#' 
#' ## The cumulative distribution and the Reliability function
#' par(mfrow=c(1, 2))
#' curve(pAPIE(x,  mu=3, sigma=0.5),
#'       from=0.0001, to=25, col="red", las=1, ylab="F(x)")
#' curve(pAPIE(x,  mu=3, sigma=0.5, lower.tail=FALSE),
#'       from=0.0001, to=25, col="red", las=1, ylab="S(x)")
#' 
#' ## The quantile function
#' p <- seq(from=0, to=0.99999, length.out=100)
#' plot(x=qAPIE(p,  mu=3, sigma=0.5), y=p, xlab="Quantile",
#'      las=1, ylab="Probability")
#' curve(pAPIE(x,  mu=3, sigma=0.5),
#'       from=0.001, add=TRUE, col="red")
#' 
#' ## The random function
#' hist(rAPIE(n=1000,  mu=3, sigma=0.5), freq=FALSE,breaks=50,
#'      xlab="x", ylim=c(0, 0.033), las=1, main="", xlim=c(0, 200))
#' curve(dAPIE(x,  mu=3, sigma=0.5),
#'       from=0.001, to=500, add=TRUE, col="red")
#' 
#' ## The Hazard function
#' curve(hAPIE(x,  mu=3, sigma=0.5), from=0.001, to=15,
#'       col="red", ylab="Hazard function", las=1)
#'
#' @references
#' \insertRef{ceren2018}{RelDists}
#'
#' @importFrom Rdpack reprompt
#'
#' @export
dAPIE <- function(x, mu, sigma, log=FALSE){
  if (any(x <= 0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(mu <= 0)) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", "")) 
  
  if (mu == 1) {
    loglik <- log(sigma) - 2 * log(x) - sigma / x
  }
  else {
    A <- log((log(mu)) / (mu - 1)) + log(sigma) - 2 * log(x) 
    B <- - (sigma / x) + exp(- sigma / x) * log(mu)
    loglik <- A + B
  }
  
  if (log == FALSE) 
    density <- exp(loglik)
  else 
    density <- loglik
  return(density)
}
#' @export
#' @rdname dAPIE
pAPIE <- function(q, mu, sigma, 
                  lower.tail=TRUE, log.p=FALSE){
  if (any(mu <= 0)) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", "")) 
  
  if (mu == 1) {
    cdf <- exp(- sigma / q)
  }
  else {
    A <-  mu^(exp(- sigma / q)) - 1
    cdf <- A / (mu - 1)
  }
  
  if (lower.tail == TRUE) 
    cdf <- cdf
  else cdf <- 1 - cdf
  if (log.p == FALSE) 
    cdf <- cdf
  else cdf <- log(cdf)
  cdf
}
#' @export
#' @rdname dAPIE
qAPIE <- function(p, mu, sigma,
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
    if (mu == 1) {
      exp(- sigma / x)
    }
    else {
      (mu^(exp(- sigma / x)) - 1) / (mu - 1)
    }
    
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
#' @rdname dAPIE
rAPIE <- function(n, mu, sigma){
  if (any(mu <= 0)) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  
  n <- ceiling(n)
  p <- runif(n)
  r <- qAPIE(p, mu, sigma)
  r
}
#' @export
#' @rdname dAPIE
hAPIE<-function(x, mu, sigma){
  if (any(x <= 0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(mu <= 0)) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))  
  
  h <- dAPIE(x, mu, sigma, log=FALSE) / 
    pAPIE(q=x, mu, sigma, lower.tail=FALSE, log.p=FALSE)
  h 
}
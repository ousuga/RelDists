#' The Exponentiated Chen distribution
#' 
#' @author Santiago Toro
#'
#' @description 
#' Density, distribution function, quantile function, 
#' random generation and hazard function for The Exponentiated Chen distribution 
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
#' The Exponentiated Chen distribution with parameters \code{mu}, 
#' \code{sigma} and \code{nu} has density given by
#' 
#' \eqn{f(x)=\mu \sigma \nu x^{\sigma - 1} \exp({x ^\sigma}) \exp\{{\nu (1 -\exp({x^\sigma})}\} 
#' \{1 - \exp({\nu (1 -\exp({x^\sigma}))})\}^{\mu-1},}
#' 
#' for \eqn{x > 0}, \eqn{\mu> 0}, \eqn{\sigma> 0}, \eqn{\nu> 0}.
#' 
#' @return 
#' \code{dEC} gives the density, \code{pEC} gives the distribution 
#' function, \code{qEC} gives the quantile function, \code{rEC}
#' generates random deviates and \code{hEC} gives the hazard function.
#'
#' @examples
#' ## The probability density function
#' curve(dEC(x, mu=35, sigma=0.25, nu=2), from=0.1, to=10,
#'       ylim=c(0, 1.5), col="red", las=1, ylab="f(x)")
#' 
#' ## The cumulative distribution and the Reliability function
#' par(mfrow=c(1, 2))
#' curve(pEC(x, mu=35, sigma=0.25, nu=2),
#'       from=0.1, to=10,  col="red", las=1, ylab="F(x)")
#' curve(pEC(x, mu=35, sigma=0.25, nu=2, lower.tail=FALSE),
#'       from=0.1, to=10, col="red", las=1, ylab="S(x)")
#' 
#' ## The quantile function
#' p <- seq(from=0, to=0.99999, length.out=100)
#' plot(x=qEC(p, mu=35, sigma=0.25, nu=2), y=p, xlab="Quantile",
#'      las=1, ylab="Probability")
#' curve(pEC(x,  mu=35, sigma=0.25, nu=2), from=0.1, add=TRUE, col="red")
#' 
#' ## The random function
#' hist(rEC(n=10000, mu=35, sigma=0.25, nu=2), freq=FALSE,
#'      xlab="x", las=1, main="")
#' curve(dEC(x, mu=35, sigma=0.25, nu=2), from=0.1, add=TRUE, col="red")
#' 
#' ## The Hazard function
#' curve(hEC(x, mu=35, sigma=0.25, nu=2), from=0.1, to=30, ylim=c(0, 1),
#'       col="red", ylab="Hazard function", las=1)
#'
#' @references
#' \insertRef{sanku2017}{RelDists}
#'
#' @importFrom Rdpack reprompt
#'
#' @export
dEC <- function(x, mu, sigma, nu, log=FALSE){
  if (any(x <= 0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(mu <= 0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu <= 0)) 
    stop(paste("nu must be positive", "\n", ""))
  
  
  term <- exp(x^sigma)
  A <- log(mu) + log(sigma) + log(nu) + (sigma - 1) * log(x) + x^sigma
  B <- nu * (1 - term) + (mu - 1) * log(1 - exp(nu * (1 - term)))
  loglik <- A + B
  
  if (log == FALSE) 
    density <- exp(loglik)
  else 
    density <- loglik
  return(density) 
}
#' @export
#' @rdname dEC
pEC <- function(q, mu, sigma, nu, 
                lower.tail=TRUE, log.p=FALSE){
  if (any(q <= 0)) 
    stop(paste("q must be positive", "\n", ""))
  if (any(mu <= 0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu <= 0)) 
    stop(paste("nu must be positive", "\n", ""))
  
  cdf <- (1 - exp(nu * (1 - exp(q^sigma))))^mu
  
  if (lower.tail == TRUE) 
    cdf <- cdf
  else cdf <- 1 - cdf
  if (log.p == FALSE) 
    cdf <- cdf
  else cdf <- log(cdf)
  cdf
}
#' @export
#' @rdname dEC
qEC <- function(p, mu, sigma, nu, 
                lower.tail=TRUE, log.p=FALSE){
  if (any(mu <= 0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu <= 0)) 
    stop(paste("nu must be positive", "\n", ""))
  
  
  if (log.p == TRUE) p <- exp(p)
  else p <- p
  if (lower.tail == TRUE) p <- p
  else p <- 1 - p
  if (any(p < 0) | any(p > 1)) 
    stop(paste("p must be between 0 and 1", "\n", ""))
  
  q <- (log(1 - (1/nu) * log(1 - p^(1/mu))))^(1/sigma)
  q
}
#' @importFrom stats runif
#' @export
#' @rdname dEC
rEC <- function(n, mu, sigma, nu){
  if(any(n <= 0))
    stop(paste("n must be positive","\n",""))
  if (any(mu <= 0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu <= 0)) 
    stop(paste("nu must be positive", "\n", ""))
  
  n <- ceiling(n)
  p <- runif(n)
  r <- qEC(p, mu, sigma, nu)
  r
}
#' @export
#' @rdname dEC
hEC <- function(x, mu, sigma, nu){
  if (any(x <= 0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(mu <= 0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu <= 0)) 
    stop(paste("nu must be positive", "\n", ""))
  
  h <- dEC(x, mu, sigma, nu, log=FALSE) / 
    pEC(x, mu, sigma, nu, lower.tail=FALSE, log.p=FALSE)
  h
}
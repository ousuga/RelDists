#' The Quasi XGamma Poisson distribution
#' 
#' @author Simon Zapata
#' 
#' @description 
#' Density, distribution function,quantile function, 
#' random generation and hazard function for the Quasi XGamma Poisson distribution 
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
#' The Quasi XGamma Poisson distribution with parameters \code{mu}, 
#' \code{sigma} and \code{nu} has density given by:
#' 
#' \eqn{f(x)= K(\mu, \sigma, \nu)(\frac {\sigma^{2} x^{2}}{2} + \mu)
#'  exp(\frac{\nu exp(-\sigma x)(1 + \mu + \sigma x + \frac {\sigma^{2}x^{2}}{2})}{1+\mu} - \sigma x),}
#' 
#' for \eqn{x > 0}, \eqn{\mu> 0}, \eqn{\sigma> 0}, \eqn{\nu> 1}.
#' 
#' where
#' 
#' \eqn{K(\mu, \sigma, \nu) = \frac{\nu \sigma}{(exp(\nu)-1)(1+\mu)}}
#' 
#' @return 
#' \code{dQXGP} gives the density, \code{pQXGP} gives the distribution 
#' function, \code{qQXGP} gives the quantile function, \code{rQXGP}
#' generates random deviates and \code{hQXGP} gives the hazard function.
#'
#' @examples  
#' ## The probability density function
#' curve(dQXGP(x, mu=0.5, sigma=1, nu=1), from=0.1, to=8,
#'       ylim=c(0, 0.6), col="red", las=1, ylab="f(x)")
#' 
#' ## The cumulative distribution and the Reliability function
#' par(mfrow=c(1, 2))
#' curve(pQXGP(x, mu=0.5, sigma=1, nu=1),
#'       from=0.1, to=8,  col="red", las=1, ylab="F(x)")
#' curve(pQXGP(x,  mu=0.5, sigma=1, nu=1, lower.tail=FALSE),
#'       from=0.1, to=8, col="red", las=1, ylab="S(x)")
#' 
#' ## The quantile function
#' p <- seq(from=0, to=0.99999, length.out=100)
#' plot(x=qQXGP(p, mu=0.5, sigma=1, nu=1), y=p, xlab="Quantile",
#'      las=1, ylab="Probability")
#' curve(pQXGP(x, mu=0.5, sigma=1, nu=1),
#'       from=0.1, add=TRUE, col="red")
#'       
#' ## The random function
#' hist(rQXGP(n=1000, mu=0.5, sigma=1, nu=1), freq=FALSE,
#'      xlab="x", ylim=c(0, 0.4), las=1, main="", xlim=c(0, 15))
#' curve(dQXGP(x, mu=0.5, sigma=1, nu=1),
#'       from=0.001, to=500, add=TRUE, col="red")
#' 
#' ## The Hazard function
#' curve(hQXGP(x, mu=0.5, sigma=1, nu=1), from=0.01, to=3,
#'       col="red", ylab="Hazard function", las=1)
#'
#' @references
#' \insertRef{subhradev2018}{RelDists}
#'
#' @importFrom Rdpack reprompt
#'
#' @export
dQXGP <- function(x, mu, sigma, nu, log = FALSE){
  if (any(x <= 0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(mu <= 0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu <= 0)) 
    stop(paste("nu must be positive", "\n", ""))
  
  K <- nu * sigma / ((exp(nu) - 1) * (1 + mu))
  A <- (1 + mu + (sigma * x) + (1 / 2 * (sigma^2) * (x^2)))
  Term <- nu * exp(-sigma * x) * A / (1+mu)
  loglik <- log(K) + log(mu + (1/2 * (sigma^2) * (x^2))) + Term  - (sigma * x)
  
  if(log == FALSE)
    density <- exp(loglik)
  else
    density <- loglik
  return(density)
}
#' @export
#' @rdname dQXGP
pQXGP <- function(q, mu, sigma, nu, lower.tail=TRUE, log.p=FALSE){
  if (any(q <= 0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(mu <= 0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu <= 0)) 
    stop(paste("nu must be positive", "\n", ""))
  
  Term1 <- nu * exp(-sigma * q)
  Term2 <- ((1 + mu + (sigma * q) + ( 1/2 * sigma^2 * q^2)) / (1 + mu))
  cdf <- (exp(nu) - exp(Term1 * Term2)) / (exp(nu) - 1)
  
  if (lower.tail == TRUE) 
    cdf <- cdf
  else cdf <- 1 - cdf
  if (log.p == FALSE) 
    cdf <- cdf
  else cdf <- log(cdf)
  cdf
}
#' @export
#' @rdname dQXGP
qQXGP <- function(p, mu, sigma, nu,
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
    
    Term1 <- nu * exp(-sigma * x)
    Term2 <- ((1 + mu + (sigma * x) + ( 1/2 * sigma^2 * x^2)) / (1 + mu))
    cdf <- (exp(nu) - exp(Term1 * Term2)) / (exp(nu) - 1)
    
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
#' @rdname dQXGP
rQXGP <- function(n, mu, sigma, nu){
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
  r <- qQXGP(p, mu, sigma, nu)
  r
}
#' @export
#' @rdname dQXGP
hQXGP <- function(x, mu, sigma, nu){
  if (any(x <= 0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(mu <= 0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu <= 0)) 
    stop(paste("nu must be positive", "\n", ""))
  
  h <- dQXGP(x, mu, sigma, nu, log=FALSE) / 
    pQXGP(q=x, mu, sigma, nu, lower.tail=FALSE, log.p=FALSE)
  h 
}

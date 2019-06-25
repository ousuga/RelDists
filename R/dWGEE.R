#' The Weighted Generalized Exponential-Exponential distribution
#' 
#' @description 
#' Density, distribution function, quantile function, 
#' random generation and hazard function for the Weighted Generalized Exponential-Exponential distribution 
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
#' The Weighted Generalized Exponential-Exponential Distribution with parameters \code{mu}, 
#' \code{sigma} and \code{nu} has density given by
#' 
#' \eqn{f(x)= \sigma \nu \exp(-\nu x) (1 - \exp(-\nu x))^{\sigma - 1} (1 - \exp(-\mu \nu x)) / 1 - \sigma B(\mu + 1, \sigma),}
#' 
#' for x > 0. 
#' 
#' @return 
#' \code{dWGEE} gives the density, \code{pWGEE} gives the distribution 
#' function, \code{qWGEE} gives the quantile function, \code{rWGEE}
#' generates random deviates and \code{hWGEE} gives the hazard function.
#'
#' @examples  
#' ## The probability density function 
#' curve(dWGEE(x, mu = 5, sigma = 0.5, nu = 1), from = 0, to = 6, 
#' ylim = c(0, 1), col = "red", las = 1, ylab = "The probability density function")
#' 
#' ## The cumulative distribution and the Reliability function
#' par(mfrow = c(1, 2))
#' curve(pWGEE(x, mu = 5, sigma = 0.5, nu = 1), from = 0, to = 6, 
#' ylim = c(0, 1), col = "red", las = 1, ylab = "The cumulative distribution function")
#' curve(pWGEE(x, mu = 5, sigma = 0.5, nu = 1, lower.tail = FALSE), 
#' from = 0, to = 6, ylim = c(0, 1), col = "red", las = 1, ylab = "The Reliability function")
#' 
#' ## The quantile function
#' p <- seq(from = 0, to = 0.99999, length.out = 100)
#' plot(x = qWGEE(p = p, mu = 5, sigma = 0.5, nu = 1), y = p, 
#' xlab = "Quantile", las = 1, ylab = "Probability")
#' curve(pWGEE(x, mu = 5, sigma = 0.5, nu = 1), from = 0, add = TRUE, 
#' col = "red")
#' 
#' ## The random function
#' hist(rWGEE(1000, mu = 5, sigma = 0.5, nu = 1), freq = FALSE, xlab = "x", 
#' ylim = c(0, 1), las = 1, main = "")
#' curve(dWGEE(x, mu = 5, sigma = 0.5, nu = 1),  from = 0, add = TRUE, 
#' col = "red", ylim = c(0, 1))
#' 
#' ## The Hazard function(
#' par(mfrow=c(1,1))
#' curve(hWGEE(x, mu = 5, sigma = 0.5, nu = 1), from = 0, to = 6, 
#' ylim = c(0, 1.4), col = "red", ylab = "The hazard function", las = 1)
#'
#' @references
#' \insertRef{Mahdavi2015}{RelDists}
#'
#' @importFrom Rdpack reprompt
#' 
#' @export
dWGEE <- function(x, mu, sigma, nu, log=FALSE) {
  if (any(x < 0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(mu <= 0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu <= 0)) 
    stop(paste("nu must be positive", "\n", ""))
  
  loglik <- log(sigma) + log(nu) - nu*x + (sigma - 1)*log(1 - exp(-nu*x)) +
    log(1 - exp(-nu*mu*x)) - log(1 - sigma*(base::beta(mu + 1, sigma)))
  
  if (log == FALSE) density <- exp(loglik)
  else density <- loglik
  return(density)
}
#' @export
#' @rdname dWGEE
pWGEE <- function(q, mu, sigma, nu, lower.tail=TRUE, log.p=FALSE){
  if (any(q < 0)) 
    stop(paste("q must be positive", "\n", ""))
  if (any(mu <= 0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu <= 0)) 
    stop(paste("nu must be positive", "\n", ""))
  
  # The incomplete beta function
  ibeta <- function(x, a, b) {
    stats::pbeta(x, a, b, lower.tail=FALSE) * base::beta(a, b)
    }
  cdf <- ((1-exp(-nu*q))^sigma - sigma * ibeta(exp(-nu*q), mu+1, sigma)) / 
    (1 - sigma * base::beta(mu + 1, sigma))
  
  if (lower.tail == TRUE) cdf <- cdf
  else cdf <- 1 - cdf
  if (log.p == FALSE) cdf <- cdf
  else cdf <- log(cdf)
  cdf 
}
#' @export
#' @rdname dWGEE
qWGEE <- function(p, mu, sigma, nu, lower.tail=TRUE, log.p=FALSE){
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
    stop(paste("p must be sigmaween 0 and 1", "\n", ""))
  
  F.inv <- function(y, mu, sigma, nu) {
    uniroot(function(x) {pWGEE(x,mu,sigma,nu) - y},
            interval=c(0, 99999))$root
  }
  F.inv <- Vectorize(F.inv)
  F.inv(p, mu, sigma, nu)
}
#' @importFrom stats runif
#' @export
#' @rdname dWGEE
rWGEE <- function(n, mu, sigma, nu){
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
  r <- qWGEE(p, mu, sigma, nu)
  r
}
#' @export
#' @rdname dWGEE
hWGEE<-function(x, mu, sigma, nu){
  if (any(x < 0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(sigma <= 0 )) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(mu <= 0)) 
    stop(paste("mu must be positive", "\n", ""))  
  if (any(nu <= 0)) 
    stop(paste("nu must be positive", "\n", ""))
  
  h <- dWGEE(x,mu, sigma,nu, log = FALSE) / 
    pWGEE(q=x,mu, sigma,nu, lower.tail=FALSE, log.p = FALSE)
  h  
}
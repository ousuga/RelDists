#' The Birnbaum-Saunders distribution
#' 
#' @description
#' Density, distribution function, quantile function, 
#' random generation and hazard function for the 
#' Birnbaum-Saunders distribution with
#' parameters \code{mu} and \code{sigma}.
#' 
#' @param x,q	vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. 
#' @param mu parameter.    
#' @param sigma parameter.
#' @param log,log.p	logical; if TRUE, probabilities p are given as log(p).	
#' @param lower.tail logical; if TRUE (default), probabilities are 
#' P[X <= x], otherwise, P[X > x].
#' 
#' @references
#' Birnbaum, Z.W. and Saunders, S.C. (1969a). A new family of life 
#' distributions. J. Appl. Prob., 6, 319-327.
#' 
#' Roquim, F. V., Ramires, T. G., Nakamura, L. R., Righetto, A. J., 
#' Lima, R. R., & Gomes, R. A. (2021). Building flexible regression 
#' models: including the Birnbaum-Saunders distribution in the 
#' gamlss package. Semina: Ciências Exatas e Tecnológicas, 
#' 42(2), 163-168.
#' 
#' @seealso \link{BS}.
#' 
#' @details 
#' The Birnbaum-Saunders with parameters \code{mu} and \code{sigma}
#' has density given by
#' 
#' \eqn{f(x) = \frac{x^{-3/2}(x+\mu)}{2\sigma\sqrt{2\pi\mu}} \exp\left(\frac{-1}{2\sigma^2}(\frac{x}{\mu}+\frac{\mu}{x}-2)\right)}
#' 
#' for \eqn{x>0}, \eqn{\mu>0} and \eqn{\sigma>0}. In this 
#' parameterization \eqn{\mu} is the median of \eqn{X}, 
#' \eqn{E(X)=\mu(1+\sigma^2/2)} and 
#' \eqn{Var(X)=(\mu\sigma)^2(1+5\sigma^2/4)}. The functions 
#' proposed here
#' corresponds to the functions created by Roquim et al. (2021)
#' with minor modifications to obtain correct log-likelihoods
#' and random samples.
#' 
#' @return 
#' \code{dBS} gives the density, \code{pBS} gives the distribution 
#' function, \code{qBS} gives the quantile function, \code{rBS}
#' generates random deviates and \code{hBS} gives the hazard function.
#' 
#' @example examples/examples_dBS.R
#' 
#' @export
dBS <- function(x, mu=1, sigma=1, log=FALSE){
  if (any(mu <= 0)) stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) stop(paste("sigma must be positive", "\n", ""))
  
  # Ensure same length vector
  ly    <- max(length(x), length(mu), length(sigma))
  xx    <- rep(x, length=ly)
  mu    <- rep(mu, length=ly)
  sigma <- rep(sigma, length=ly)
  
  # Temporal change for invalid x's
  xx[x <= 0] <- 0.5
  xx[is.infinite(x)] <- 0.5
  
  # pdf in log-scale
  p <- -1.5*log(xx)+log(xx+mu)-log(2*sigma)-0.5*log(2*pi*mu)-0.5*(xx/mu+mu/xx-2)/sigma^2
  
  # Assign values for invalid x's
  p[x <= 0] <- -Inf
  p[is.infinite(x)] <- -Inf
  
  if (log == FALSE)
    p <- exp(p)
  
  return(p)
}
#' @export
#' @importFrom stats pnorm
#' @rdname dBS
pBS <- function(q, mu=1, sigma=1, lower.tail=TRUE, log.p=FALSE){
  
  if (any(mu <= 0))    stop("parameter mu has to be positive!")
  if (any(sigma <= 0)) stop("parameter sigma has to be positive!")
  
  # Ensure same length vector
  ly    <- max(length(q), length(mu), length(sigma))
  qq    <- rep(q, length=ly)
  mu    <- rep(mu, length=ly)
  sigma <- rep(sigma, length=ly)
  
  # Temporal change for invalid x's
  qq[q <= 0] <- 0.5
  qq[q == Inf] <- 0.5
  
  # The cumulative
  cdf <- pnorm(((qq/mu)^0.5-(mu/qq)^0.5)/sigma)
  
  # Assign values for invalid x's
  cdf[q <= 0] <- 0
  cdf[q == Inf] <- 1
  
  if (lower.tail == FALSE)
    cdf <- 1 - cdf
  if (log.p == TRUE)
    cdf <- log(cdf)
  
  return(cdf)
}
#' @importFrom stats uniroot qnorm
#' @export
#' @rdname dBS
qBS <- function(p, mu=1, sigma=1, lower.tail = TRUE, log.p = FALSE){
  if (any(mu <= 0)) stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) stop(paste("sigma must be positive", "\n", ""))

  # To adjust the probability
  if (log.p == TRUE)
    p <- exp(p)
  if (lower.tail == FALSE)
    p <- 1 - p
  
  # Ensure same length vector
  ly <- max(length(p), length(mu), length(sigma))
  pp <- rep(p, length=ly)
  mu <- rep(mu, length=ly)
  sigma <- rep(sigma, length=ly)
  
  # Temporal change for invalid p's
  pp[p < 0]  <-  0.5
  pp[p > 1]  <-  0.5
  pp[p == 1] <-  0.5
  pp[p == 0] <-  0.5
  
  # The quantile
  w <- sigma * qnorm(pp)/2
  q <- mu * (w + sqrt(w^2+1))^2
  
  # To deal with invalid p's
  q[p <  0] <- NaN
  q[p >  1] <- NaN
  q[p == 1] <- Inf
  q[p == 0] <- 0
  
  return(q)
}
#' @importFrom stats runif
#' @export
#' @rdname dBS
rBS <- function(n, mu=1, sigma=1){
  if (any(mu <= 0))     stop("parameter mu has to be positive!")
  if (any(sigma <= 0))  stop("parameter sigma has to be positive!")
  if (any(n <= 0))      stop(paste("n must be a positive integer", "\n", ""))
  
  n <- ceiling(n)
  u <- runif(n=n)
  x <- qBS(p=u, mu=mu, sigma=sigma)
  return(x)
}
#' @export
#' @rdname dBS
hBS <- function(x, mu, sigma){
  if (any(x < 0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(mu <= 0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  
  h <- dBS(x, mu, sigma) / pBS(x, mu, sigma, lower.tail=FALSE)
  h
}
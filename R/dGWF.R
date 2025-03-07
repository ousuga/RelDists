#' Generalized Weibull distribution
#'
#' @author Jaime Mosquera, \email{jmosquerag@unal.edu.co}
#'
#' @description
#' Density, distribution function, quantile function, random generation and
#' hazard function for the generalized Weibull distribution with parameters
#' \code{mu}, \code{sigma} and \code{nu}.
#'
#' @param x,q	vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations.
#' @param mu parameter one.
#' @param sigma parameter two.
#' @param nu parameter three.
#' @param log,log.p	logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE (default), probabilities are
#' P[T <= t], otherwise, P[T > t].
#'
#' @references
#' Mudholkar, G. S., & Kollia, G. D. (1994). Generalized Weibull family: a structural analysis. Communications in statistics-theory and methods, 23(4), 1149-1171.
#'
#' @details
#' The generalized Weibull with parameters \code{mu}, \code{sigma} and
#' \code{nu} has density given by
#'
#' \deqn{f(x) = \mu \sigma x^{\sigma - 1}
#'    \left( 1 - \mu \nu x^\sigma \right)^{\frac{1}{\nu} - 1}}
#'
#' for \eqn{x > 0}, \eqn{\mu>0}, \eqn{\sigma>0} and \eqn{-\infty<\nu<\infty}.
#'
#' @return
#' \code{dGWF} gives the density, \code{pGWF} gives the distribution
#' function, \code{qGWF} gives the quantile function, \code{rGWF}
#' generates random deviates and \code{hGWF} gives the hazard function.
#'
#' @example examples/examples_dGWF.R
#'
#' @export
dGWF <- function(x, mu, sigma, nu, log = FALSE) {
  if (any(x < 0)) {
    stop(paste("x must be positive", "\n", ""))
  }
  if (mu <= 0 || sigma <= 0) {
    stop(paste("mu, sigma, and nu must be positive", "\n", ""))
  }
  f <- mu * sigma * x^(sigma - 1) * (1 - mu * nu * x^sigma)^(1 / nu - 1)
  if (log == FALSE) {
    dens <- f
  } else {
    dens <- log(f)
  }
  return(dens)
}
#' @export
#' @rdname dGWF
pGWF <- function(q, mu, sigma, nu, lower.tail = TRUE, log.p = FALSE) {
  if (any(q < 0)) {
    stop(paste("q must be positive", "\n", ""))
  }
  if (mu <= 0 || sigma <= 0) {
    stop(paste("mu, sigma, and nu must be positive", "\n", ""))
  }
  cdf <- 1 - (1 - mu * nu * q^sigma)^(1 / nu)
  if (lower.tail == TRUE) {
    cdf <- cdf
  } else {
    cdf <- 1 - cdf
  }
  if (log.p == FALSE) {
    cdf <- cdf
  } else {
    cdf <- log(cdf)
  }
  return(cdf)
}
#' @export
#' @rdname dGWF
qGWF <- function(p, mu, sigma, nu) {
  if (mu <= 0 || sigma <= 0) {
    stop(paste("mu, sigma, and nu must be positive", "\n", ""))
  }
  if (any(p < 0) | any(p > 1)) {
    stop(paste("p must be between 0 and 1", "\n", ""))
  }
  q <- ((1 - (1 - p)^nu) / (mu * nu))^(1 / sigma)
  return(q)
}
#' @export
#' @rdname dGWF
rGWF <- function(n, mu, sigma, nu) {
  if (mu <= 0 || sigma <= 0) {
    stop(paste("mu, sigma, and nu must be positive", "\n", ""))
  }
  if (n <= 0) {
    stop(paste("n must be positive", "\n", ""))
  }
  p <- runif(n)
  q <- qGWF(p, mu, sigma, nu)
  return(q)
}
#' @export
#' @rdname dGWF
hGWF <- function(x, mu, sigma, nu) {
  if (any(x < 0)) {
    stop(paste("x must be positive", "\n", ""))
  }
  if (mu <= 0 || sigma <= 0) {
    stop(paste("mu, sigma, and nu must be positive", "\n", ""))
  }
  h <- mu * sigma * x^(sigma - 1) * (1 - mu * nu * x^sigma)^(-1)
  return(h)
}

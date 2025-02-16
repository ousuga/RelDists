#' The New Modified Weibull distribution
#' @description
#' Density, distribution function, quantile function,
#' random generation and hazard function for the reduced new modified Weibull
#' distribution with parameters \code{mu}, \code{sigma} and \code{nu}.
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
#' @details
#' The reduced new modified Weibull with parameters \code{mu}, \code{sigma}
#' and \code{theta} has density given by
#'
#' \deqn{f(x) = \frac{1}{2 \sqrt{x}}
#'    \left( \mu + \sigma (1 + 2 \nu x) e^{\nu x} \right)
#'    e^{-\mu \sqrt{x} - \sigma \sqrt{x} e^{\nu x}}}
#'
#' for x > 0, \deqn{\mu}, \deqn{\sigma} and \deqn{\nu} > 0.
#'
#' @return
#' \code{dRNMW} gives the density, \code{pRNMW} gives the distribution
#' function, \code{qRNMW} gives the quantile function, \code{rRNMW}
#' generates random deviates and \code{hRNMW} gives the hazard function.
#'
#' @example examples/examples_dRNMW.R
#'
#' @references
#' \insertRef{almalki2013reduced}{RelDists}
#'
#' @export
dRNMW <- function(x, mu, sigma, nu, log = FALSE) {
  if (any(x < 0)) {
    stop(paste("x must be positive", "\n", ""))
  }
  if (mu <= 0 || sigma <= 0 || nu <= 0) {
    stop(paste("mu, sigma, and nu must be positive", "\n", ""))
  }
  e_term1 <- exp(nu * x)
  e_term2 <- exp(-mu * sqrt(x) - sigma * sqrt(x) * exp(nu * x))
  f <- (1 / (2 * sqrt(x))) * (mu + sigma * (1 + 2 * nu * x) * e_term1) * e_term2
  if (log == FALSE) {
    dens <- f
  } else {
    dens <- log(f)
  }
  return(dens)
}
#' @export
#' @rdname dRNMW
pRNMW <- function(
    q, mu, sigma, nu, lower.tail = TRUE, log.p = FALSE) {
  if (any(q < 0)) {
    stop(paste("q must be positive", "\n", ""))
  }
  if (mu <= 0 || sigma <= 0 || nu <= 0) {
    stop(paste("mu, sigma, and nu must be positive", "\n", ""))
  }
  cdf <- 1 - exp(-mu * sqrt(q) - sigma * sqrt(q) * exp(nu * q))
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
#' @rdname dRNMW
qRNMW <- function(p, mu, sigma, nu, lower.tail = TRUE, log.p = FALSE) {
  if (mu <= 0 || sigma <= 0 || nu <= 0) {
    stop(paste("mu, sigma, and nu must be positive", "\n", ""))
  }
  if (any(p < 0) | any(p > 1)) {
    stop(paste("p must be between 0 and 1", "\n", ""))
  }

  if (lower.tail == TRUE) {
    p <- p
  } else {
    p <- 1 - p
  }
  if (log.p == FALSE) {
    p <- p
  } else {
    p <- exp(p)
  }
  fda <- function(x, mu, sigma, nu) {
    1 - exp(-mu * sqrt(x) - sigma * sqrt(x) * exp(nu * x))
  }

  fda1 <- function(x, mu, sigma, nu, p) {
    fda(x, mu, sigma, nu) - p
  }

  r_de_la_funcion <- function(mu, sigma, nu, p) {
    uniroot(fda1, interval = c(0, 1e+06), mu, sigma, nu, p)$root
  }

  r_de_la_funcion <- Vectorize(r_de_la_funcion)
  q <- r_de_la_funcion(mu, sigma, nu, p)
  return(q)
}
#' @export
#' @rdname dRNMW
rRNMW <- function(n, mu, sigma, nu) {
  if (any(n <= 0)) {
    stop(paste("n must be positive", "\n", ""))
  }
  if (mu <= 0 || sigma <= 0 || nu <= 0) {
    stop(paste("mu, sigma, and nu must be positive", "\n", ""))
  }
  p <- runif(n)
  r <- qRNMW(p, mu, sigma, nu)
  return(r)
}
#' @export
hRNMW <- function(x, mu, sigma, nu) {
  if (any(x < 0)) {
    stop(paste("x must be positive", "\n", ""))
  }
  if (mu <= 0 || sigma <= 0 || nu <= 0) {
    stop(paste("mu, sigma, and nu must be positive", "\n", ""))
  }
  f <- dRNMW(x, mu, sigma, nu, log = FALSE)
  S <- pRNMW(
    q = x, mu, sigma, nu, lower.tail = FALSE, log.p = FALSE
  )
  h <- f / S
  return(h)
}

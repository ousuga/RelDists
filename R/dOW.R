#' The Odd Weibull Distribution
#'
#' @author Jaime Mosquera Gutiérrez \email{jmosquerag@unal.edu.co}
#'
#' @description
#' Density, distribution function, quantile function,
#' random generation and hazard function for the Odd Weibull distribution with
#' parameters \code{mu}, \code{sigma} and \code{nu}.
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
#' The Odd Weibull with parameters \code{mu}, \code{sigma}
#' and \code{nu} has density given by
#'
#' \eqn{f(x) = \left( \frac{\sigma\nu}{x} \right) (\mu x)^\sigma
#'      e^{(\mu x)^\sigma} \left(e^{(\mu x)^{\sigma}}-1\right)^{\nu-1}
#'      \left[ 1 + \left(e^{(\mu x)^{\sigma}}-1\right)^\nu \right]^{-2}}
#'
#' for x > 0.
#'
#' @return
#' \code{dOW} gives the density, \code{pOW} gives the distribution
#' function, \code{qOW} gives the quantile function, \code{rOW}
#' generates random deviates and \code{hOW} gives the hazard function.
#'
#' @example examples/examples_dOW.R
#'
#' @references
#' Cooray, K. (2006). Generalization of the Weibull distribution: 
#' the odd Weibull family. Statistical Modelling, 6(3), 265-277.
#'
#' @export
dOW <- function(x, mu, sigma, nu, log = FALSE) {
  if (any(x < 0)) {
    stop(paste("x must be positive", "\n", ""))
  }
  if (any(mu <= 0)) {
    stop(paste("mu must be positive", "\n", ""))
  }
  if (any(sigma * nu <= 0)) {
    stop(paste("Product sigma*nu must be positive", "\n", ""))
  }

  prod1 <- (mu * x)^sigma
  loglik <- log(sigma * nu) - log(x) + sigma * (log(mu) + log(x)) +
    (mu * x)^sigma + (nu - 1) * log(expm1(prod1)) -
    2 * log1p((expm1(prod1))^nu)

  if (log == FALSE) {
    dens <- exp(loglik)
  } else {
    dens <- loglik
  }
  return(dens)
}

#' @export
#' @rdname dOW
pOW <- function(q, mu, sigma, nu, lower.tail = TRUE, log.p = FALSE) {
  if (any(q < 0)) {
    stop(paste("q must be positive", "\n", ""))
  }
  if (any(mu <= 0)) {
    stop(paste("mu must be positive", "\n", ""))
  }
  if (any(sigma * nu <= 0)) {
    stop(paste("Product sigma*nu must be positive", "\n", ""))
  }

  prod1 <- (mu * q)^sigma
  cdf <- 1 - (1 + (expm1(prod1))^nu)^(-1)

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
#' @rdname dOW
qOW <- function(p, mu, sigma, nu, lower.tail = TRUE, log.p = FALSE) {
  if (any(mu <= 0)) {
    stop(paste("mu must be positive", "\n", ""))
  }
  if (any(sigma * nu <= 0)) {
    stop(paste("Product sigma*nu must be positive", "\n", ""))
  }
  if (any(p < 0) | any(p > 1)) 
    stop(paste("p must be between 0 and 1", "\n", ""))
  
  if (log.p == TRUE) {
    p <- exp(p)
  } else {
    p <- p
  }
  if (lower.tail == TRUE) {
    p <- p
  } else {
    p <- 1 - p
  }
  q <- (1 / mu) * (log1p((p * (1 - p)^(-1))^(1 / nu)))^(1 / sigma)
  return(q)
}

#' @export
#' @rdname dOW
rOW <- function(n, mu, sigma, nu) {
  if (any(n <= 0)) {
    stop(paste("n must be positive", "\n", ""))
  }
  if (any(mu <= 0)) {
    stop(paste("mu must be positive", "\n", ""))
  }
  if (any(sigma * nu <= 0)) {
    stop(paste("Product sigma*nu must be positive", "\n", ""))
  }

  n <- ceiling(n)
  p <- runif(n)
  r <- qOW(p, mu, sigma, nu)
  return(r)
}
#' @export
#' @rdname dOW
hOW <- function(x, mu, sigma, nu) {
  if (any(x < 0)) {
    stop(paste("x must be positive", "\n", ""))
  }
  if (any(mu <= 0)) {
    stop(paste("mu must be positive", "\n", ""))
  }
  if (any(sigma * nu <= 0)) {
    stop(paste("Product sigma*nu must be positive", "\n", ""))
  }

  h <- dOW(x, mu, sigma, nu, log = FALSE) / pOW(
    q = x, mu, sigma, nu,
    lower.tail = FALSE, log.p = FALSE
  )
  return(h)
}

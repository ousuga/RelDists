#' Lindley distribution
#' 
#' @author Freddy Hernandez, \email{fhernanb@@unal.edu.co}
#' 
#' @description 
#' Density, distribution function, quantile function, 
#' random generation and hazard function for the Lindley distribution
#' with parameter \code{mu}.
#' 
#' @param x,q	vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. 
#' @param mu parameter.
#' @param log,log.p	logical; if TRUE, probabilities p are given as log(p).	
#' @param lower.tail logical; if TRUE (default), probabilities are P[X <= x], otherwise, P[X > x].
#' 
#' @details 
#' Lindley Distribution with parameter \code{mu} has density given by
#' 
#' \eqn{f(x) = \frac{\mu^2}{\mu+1} (1+x) \exp(-\mu x),}
#' 
#' for x > 0 and \eqn{\mu > 0}. These function were taken form LindleyR package.
#' 
#' @return 
#' \code{dLIN} gives the density, \code{pLIN} gives the distribution 
#' function, \code{qLIN} gives the quantile function, \code{rLIN}
#' generates random deviates and \code{hLIN} gives the hazard function.
#'
#' @example examples/examples_dLIN.R  
#'
#' @references
#' Lindley, D. V. (1958). Fiducial distributions and Bayes' theorem. 
#' Journal of the Royal Statistical Society. 
#' Series B (Methodological), 102-107.
#'
#' @export
dLIN <- function (x, mu, log = FALSE) {
  stopifnot(mu > 0)
  if (log) {
    t1 <- log(mu)
    t4 <- log1p(mu)
    t6 <- log1p(x)
    -mu * x + 2 * t1 - t4 + t6
  }
  else {
    t1 <- mu^2
    t7 <- exp(-mu * x)
    t1/(1 + mu) * (1 + x) * t7
  }
}
#' @export
#' @rdname dLIN
pLIN <- function (q, mu, lower.tail = TRUE, log.p = FALSE) {
  stopifnot(mu > 0)
  if (lower.tail) {
    t1 <- mu * q
    t6 <- exp(-t1)
    cdf <- 1 - (1 + t1/(1 + mu)) * t6
  }
  else {
    t1 <- mu * q
    t6 <- exp(-t1)
    cdf <- (1 + t1/(1 + mu)) * t6
  }
  if (log.p) 
    return(log(cdf))
  else return(cdf)
}
#' @importFrom lamW lambertWm1
#' @export
#' @rdname dLIN
qLIN <- function (p, mu, lower.tail = TRUE, log.p = FALSE) 
{
  stopifnot(mu > 0)
  if (lower.tail) {
    t1 <- 1 + mu
    t4 <- exp(-t1)
    t6 <- lambertWm1(t1 * (p - 1) * t4)
    qtf <- -(t6 + 1 + mu)/mu
  }
  else {
    t1 <- 1 + mu
    t3 <- exp(-t1)
    t5 <- lambertWm1(-p * t1 * t3)
    qtf <- -(t5 + 1 + mu)/mu
  }
  if (log.p) 
    return(log(qtf))
  else return(qtf)
}
#' @importFrom stats runif rbinom rgamma
#' @export
#' @rdname dLIN
rLIN <- function (n, mu) 
{
  stopifnot(mu > 0)
  x <- rbinom(n, size = 1, prob = mu/(1 + mu))
  x <- x * rgamma(n, shape = 1, rate = mu) + (1 - x) * rgamma(n, shape = 2, rate = mu)
  x
}
#' @export
#' @rdname dLIN
hLIN <- function(x, mu, log = FALSE)
{
  stopifnot(mu > 0)
  if(log)
  {
    t1 <- log(mu)
    t5 <- log1p(mu * x + mu)
    t7 <- log1p(x)
    2 * t1 - t5 + t7
  }
  else
  {
    t1 <- mu ^ 2
    t1 / (mu * x + mu + 1) * (1 + x)
  }
}


#' The Lindley family
#' 
#' @author Freddy Hernandez \email{fhernanb@unal.edu.co}
#' 
#' @description 
#' The function \code{LIN()} defines the Lindley distribution with only one parameter 
#' for a \code{gamlss.family} object to be used in GAMLSS fitting 
#' using the function \code{gamlss()}.
#' 
#' @param mu.link defines the mu.link, with "log" link as the default for the mu parameter.
#' 
#' @details 
#' The Lindley with parameter \code{mu} has density given by
#' 
#' \eqn{f(x) = \frac{\mu^2}{\mu+1} (1+x) \exp(-\mu x),}
#'
#' for x > 0 and \eqn{\mu > 0}.
#' 
#' @returns Returns a gamlss.family object which can be used to fit a LIN distribution in the \code{gamlss()} function.
#' 
#' @example examples/examples_LIN.R
#' 
#' @references
#' Lindley, D. V. (1958). Fiducial distributions and Bayes' theorem. 
#' Journal of the Royal Statistical Society. 
#' Series B (Methodological), 102-107.
#' 
#' @importFrom gamlss.dist checklink
#' @importFrom gamlss rqres.plot
#' @export
LIN <- function (mu.link = "log") {
  mstats <- checklink("mu.link", "Lindely", substitute(mu.link), 
                      c("inverse", "log", "sqrt", "identity"))
  structure(list(family = c("LIN", "Lindley"), parameters = list(mu = TRUE), 
                 nopar = 1, 
                 type = "Continuous", 
                 mu.link = as.character(substitute(mu.link)), 
                 mu.linkfun = mstats$linkfun, 
                 mu.linkinv = mstats$linkinv, 
                 mu.dr = mstats$mu.eta, 
                 dldm = function(y, mu) 2 / mu - 1 / (mu + 1) - y, 
                 d2ldm2 = function(mu) 1 / (mu + 1)^2 - 2 / mu^2, 
                 G.dev.incr = function(y, mu, ...) -2 * dLIN(x = y, mu = mu, log = TRUE), 
                 rqres = expression(rqres(pfun = "pLIN", 
                                          type = "Continuous",
                                          y = y, 
                                          mu = mu)), 
                 mu.initial = expression(mu <- rep((-(mean(y)-1) + sqrt((mean(y)-1)^2 + 8 * mean(y))) / (2 * mean(y)), length(y))),
                 mu.valid = function(mu) all(mu > 0), 
                 y.valid = function(y) all(y > 0), 
                 mean = function(mu) (mu + 2) / (mu * (mu + 1)), 
                 variance = function(mu) 2 * (mu + 3) / (mu^2 * (mu + 1)) - ((mu + 2) / (mu * (mu + 1)))^2
  ), 
  class = c("gamlss.family", "family"))
}


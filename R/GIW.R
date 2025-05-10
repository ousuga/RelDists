#' The Generalized Inverse Weibull family
#' 
#' @author Amylkar Urrea Montoya, \email{amylkar.urrea@@udea.edu.co}
#' 
#' @description 
#' The Generalized Inverse Weibull family
#' 
#' @param mu.link defines the mu.link, with "log" link as the default for the mu parameter.
#' @param sigma.link defines the sigma.link, with "log" link as the default for the sigma.
#' @param nu.link defines the nu.link, with "log" link as the default for the nu parameter.
#' 
#' @seealso \link{dGIW}
#' 
#' @details 
#' The Generalized Inverse Weibull distribution with parameters \code{mu}, 
#' \code{sigma} and \code{nu} has density given by
#' 
#' \eqn{f(x) = \nu \sigma \mu^{\sigma} x^{-(\sigma + 1)} exp \{-\nu (\frac{\mu}{x})^{\sigma}\},}
#' 
#' for x > 0. 
#' 
#' @returns Returns a gamlss.family object which can be used to fit a GIW distribution in the \code{gamlss()} function.
#' 
#' @example examples/examples_GIW.R 
#' 
#' @references
#' Almalki, S. J., & Nadarajah, S. (2014). Modifications of the 
#' Weibull distribution: A review. Reliability Engineering & 
#' System Safety, 124, 32-55.
#' 
#' De Gusmao, F. R., Ortega, E. M., & Cordeiro, G. M. (2011). 
#' The generalized inverse Weibull distribution. Statistical 
#' Papers, 52, 591-619.
#'
#' @importFrom gamlss.dist checklink
#' @importFrom gamlss rqres.plot
#' @export
GIW <- function (mu.link="log", sigma.link="log", nu.link="log") {
  mstats <- checklink("mu.link", "Generalized Inverse Weibull", 
                      substitute(mu.link), c("log", "own"))
  dstats <- checklink("sigma.link", "Generalized Inverse Weibull",
                      substitute(sigma.link), c("log", "own"))
  vstats <- checklink("nu.link", "Generalized Inverse Weibull", 
                      substitute(nu.link), c("log", "own"))
  
  structure(list(family=c("GIW", "Generalized Inverse Weibull"), 
                 parameters=list(mu=TRUE, sigma=TRUE, nu=TRUE), 
                 nopar=3, 
                 type="Continuous", 
                 
                 mu.link = as.character(substitute(mu.link)), 
                 sigma.link = as.character(substitute(sigma.link)), 
                 nu.link = as.character(substitute(nu.link)),
                 
                 mu.linkfun = mstats$linkfun, 
                 sigma.linkfun = dstats$linkfun, 
                 nu.linkfun = vstats$linkfun,
                 
                 mu.linkinv = mstats$linkinv, 
                 sigma.linkinv = dstats$linkinv, 
                 nu.linkinv = vstats$linkinv,
                 
                 mu.dr = mstats$mu.eta, 
                 sigma.dr = dstats$mu.eta, 
                 nu.dr = vstats$mu.eta,
                 
                 # Primeras derivadas ---------------------------------
                 dldm = function(y, mu, sigma, nu) {
                   A    <- (nu * sigma * (mu / y)^(sigma - 1)) / y
                   dldm <-  sigma / mu - A
                   dldm
                 },
                 
                 dldd = function(y, mu, sigma, nu) {
                   B    <- nu * (mu / y)^sigma * log(mu / y) 
                   dldd <- 1 / sigma + log(mu) - log(y) - B
                   dldd
                 },
                 
                 dldv = function(y, mu, sigma, nu) {
                   dldv <- 1 / nu - (mu / y)^sigma
                   dldv
                 },
                 
                 
                 # Segundas derivadas ---------------------------------
                 d2ldm2 = function(y, mu, sigma, nu) {
                   A    <- (nu * sigma * (mu / y)^(sigma - 1)) / y
                   dldm <-  sigma / mu - A
                   d2ldm2 <- -dldm * dldm
                   d2ldm2
                 },
                 
                 d2ldmdd = function(y, mu, sigma, nu) {
                   A    <- (nu * sigma * (mu / y)^(sigma - 1)) / y
                   dldm <-  sigma / mu - A
                   B    <- nu * (mu / y)^sigma * log(mu / y) 
                   dldd <- 1 / sigma + log(mu) - log(y) - B
                   d2ldmdd <- -dldm * dldd
                   d2ldmdd
                 },
                 
                 d2ldmdv = function(y, mu, sigma, nu) {
                   A    <- (nu * sigma * (mu / y)^(sigma - 1)) / y
                   dldm <-  sigma / mu - A
                   dldv <- 1 / nu - (mu / y)^sigma
                   d2ldmdv <- -dldm * dldv
                   d2ldmdv
                 },
                 
                 d2ldd2 = function(y, mu, sigma, nu) {
                   B    <- nu * (mu / y)^sigma * log(mu / y) 
                   dldd <- 1 / sigma + log(mu) - log(y) - B
                   d2ldd2 <- -dldd * dldd
                   d2ldd2
                 },
                 
                 d2ldddv = function(y, mu, sigma, nu) {
                   B    <- nu * (mu / y)^sigma * log(mu / y) 
                   dldd <- 1 / sigma + log(mu) - log(y) - B
                   dldv <- 1 / nu - (mu / y)^sigma
                   d2ldddv <- -dldd * dldv
                   d2ldddv
                 },
                 
                 d2ldv2 = function(y, mu, sigma, nu) {
                   dldv <- 1 / nu - (mu / y)^sigma
                   d2ldv2 <- -dldv * dldv
                   d2ldv2
                 },
                 
                 G.dev.incr = function(y, mu, sigma, nu, ...) -2*dGIW(y, mu, sigma, nu, log=TRUE), 
                 rqres = expression(rqres(pfun="pGIW", type="Continuous", y=y, mu=mu, sigma=sigma, nu=nu)), 
                 
                 mu.initial = expression(mu       <- rep(1, length(y))), 
                 sigma.initial = expression(sigma <- rep(1, length(y))), 
                 nu.initial = expression(nu       <- rep(1, length(y))),
                 
                 mu.valid = function(mu)       all(mu > 0), 
                 sigma.valid = function(sigma) all(sigma > 0), 
                 nu.valid = function(nu)       all(nu > 0),
                 
                 y.valid = function(y) all(y > 0)
  ), 
  class=c("gamlss.family", "family"))
}

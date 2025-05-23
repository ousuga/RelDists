#' The Marshall-Olkin Extended Weibull family
#' 
#' @author Amylkar Urrea Montoya, \email{amylkar.urrea@@udea.edu.co}
#' 
#' @description 
#' The Marshall-Olkin Extended Weibull family
#' 
#' @param mu.link defines the mu.link, with "log" link as the default for the mu parameter.
#' @param sigma.link defines the sigma.link, with "log" link as the default for the sigma.
#' @param nu.link defines the nu.link, with "log" link as the default for the nu parameter.
#' 
#' @seealso \link{dMOEW}
#' 
#' @details 
#' The Marshall-Olkin Extended Weibull distribution with parameters \code{mu}, 
#' \code{sigma} and \code{nu} has density given by
#' 
#' \eqn{f(x) = \frac{\mu \sigma \nu (\nu x)^{\sigma - 1} exp\{{-(\nu x )^{\sigma}}\}}{\{1-(1-\mu) exp\{{-(\nu x )^{\sigma}}\} \}^{2}},}
#' 
#' for x > 0. 
#' 
#' @returns Returns a gamlss.family object which can be used to fit a MOEW distribution in the \code{gamlss()} function.
#' 
#' @example examples/examples_MOEW.R 
#' 
#' @references
#' Almalki, S. J., & Nadarajah, S. (2014). Modifications of the 
#' Weibull distribution: A review. Reliability Engineering & 
#' System Safety, 124, 32-55.
#' 
#' Ghitany, M. E., Al-Hussaini, E. K., & Al-Jarallah, R. A. (2005). 
#' Marshall–Olkin extended Weibull distribution and its application 
#' to censored data. Journal of Applied Statistics, 32(10), 1025-1034.
#'
#' @importFrom gamlss.dist checklink
#' @importFrom gamlss rqres.plot
#' @export
MOEW <- function (mu.link="log", sigma.link="log", nu.link="log") {
  mstats <- checklink("mu.link", "Marshall-Olkin Extended Weibull", 
                      substitute(mu.link), c("log", "own"))
  dstats <- checklink("sigma.link", "Marshall-Olkin Extended Weibull",
                      substitute(sigma.link), c("log", "own"))
  vstats <- checklink("nu.link", "Marshall-Olkin Extended Weibull", 
                      substitute(nu.link), c("log", "own"))
  
  structure(list(family=c("MOEW", "Marshall-Olkin Extended Weibull"), 
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
                   A    <- exp(-(nu * y)^sigma)
                   B    <- 1 - (1 - mu) * exp(-(nu * y)^sigma)
                   dldm <-  1 / mu - 2 * A / B
                   dldm
                 },
                 
                 dldd = function(y, mu, sigma, nu) {
                   C    <- (1 / sigma) - log(nu * y) * (nu * y)^sigma + log(nu*y)
                   D    <- (1 - mu) * exp(-(nu * y)^sigma) * log(nu * y) * (nu * y)^sigma
                   E    <- 1 - (1 - mu) * exp(-(nu * y)^sigma)
                   dldd <- C + (-2 * D / E)
                   dldd
                 },
                 
                 dldv = function(y, mu, sigma, nu) {
                   G    <- 1 / nu + (sigma - 1) / nu - sigma * y * (nu * y)^(sigma - 1)
                   H    <- - 2 * (1 - mu) * sigma * exp(-(nu * y)^sigma) * y * (nu * y)^(sigma - 1)
                   I    <- 1 - (1 - mu) * exp(-(nu * y)^sigma)
                   dldv <- G + H / I
                   dldv
                 },
                 
                 
                 # Segundas derivadas ---------------------------------
                 d2ldm2 = function(y, mu, sigma, nu) {
                   A    <- exp(-(nu * y)^sigma)
                   B    <- 1 - (1 - mu) * exp(-(nu * y)^sigma)
                   dldm <-  1 / mu - 2 * A / B
                   d2ldm2 <- -dldm * dldm
                   d2ldm2
                 },
                 
                 d2ldmdd = function(y, mu, sigma, nu) {
                   A    <- exp(-(nu * y)^sigma)
                   B    <- 1 - (1 - mu) * exp(-(nu * y)^sigma)
                   dldm <-  1 / mu - 2 * A / B
                   C    <- (1 / sigma) - log(nu * y) * (nu * y)^sigma + log(nu*y)
                   D    <- (1 - mu) * exp(-(nu * y)^sigma) * log(nu * y) * (nu * y)^sigma
                   E    <- 1 - (1 - mu) * exp(-(nu * y)^sigma)
                   dldd <- C + (-2 * D / E)
                   d2ldmdd <- -dldm * dldd
                   d2ldmdd
                 },
                 
                 d2ldmdv = function(y, mu, sigma, nu) {
                   A    <- exp(-(nu * y)^sigma)
                   B    <- 1 - (1 - mu) * exp(-(nu * y)^sigma)
                   dldm <-  1 / mu - 2 * A / B
                   G    <- 1 / nu + (sigma - 1) / nu - sigma * y * (nu * y)^(sigma - 1)
                   H    <- - 2 * (1 - mu) * sigma * exp(-(nu * y)^sigma) * y * (nu * y)^(sigma - 1)
                   I    <- 1 - (1 - mu) * exp(-(nu * y)^sigma)
                   dldv <- G + H / I
                   d2ldmdv <- -dldm * dldv
                   d2ldmdv
                 },
                 
                 d2ldd2 = function(y, mu, sigma, nu) {
                   C    <- (1 / sigma) - log(nu * y) * (nu * y)^sigma + log(nu*y)
                   D    <- (1 - mu) * exp(-(nu * y)^sigma) * log(nu * y) * (nu * y)^sigma
                   E    <- 1 - (1 - mu) * exp(-(nu * y)^sigma)
                   dldd <- C + (-2 * D / E)
                   d2ldd2 <- -dldd * dldd
                   d2ldd2
                 },
                 
                 d2ldddv = function(y, mu, sigma, nu) {
                   C    <- (1 / sigma) - log(nu * y) * (nu * y)^sigma + log(nu*y)
                   D    <- (1 - mu) * exp(-(nu * y)^sigma) * log(nu * y) * (nu * y)^sigma
                   E    <- 1 - (1 - mu) * exp(-(nu * y)^sigma)
                   dldd <- C + (-2 * D / E)
                   G    <- 1 / nu + (sigma - 1) / nu - sigma * y * (nu * y)^(sigma - 1)
                   H    <- - 2 * (1 - mu) * sigma * exp(-(nu * y)^sigma) * y * (nu * y)^(sigma - 1)
                   I    <- 1 - (1 - mu) * exp(-(nu * y)^sigma)
                   dldv <- G + H / I
                   d2ldddv <- -dldd * dldv
                   d2ldddv
                 },
                 
                 d2ldv2 = function(y, mu, sigma, nu) {
                   G    <- 1 / nu + (sigma - 1) / nu - sigma * y * (nu * y)^(sigma - 1)
                   H    <- - 2 * (1 - mu) * sigma * exp(-(nu * y)^sigma) * y * (nu * y)^(sigma - 1)
                   I    <- 1 - (1 - mu) * exp(-(nu * y)^sigma)
                   dldv <- G + H / I
                   d2ldv2 <- -dldv * dldv
                   d2ldv2
                 },
                 
                 G.dev.incr = function(y, mu, sigma, nu, ...) -2*dMOEW(y, mu, sigma, nu, log=TRUE), 
                 rqres = expression(rqres(pfun="pMOEW", type="Continuous", y=y, mu=mu, sigma=sigma, nu=nu)), 
                 
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

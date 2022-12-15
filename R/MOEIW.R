#' The Marshall-Olkin Extended Inverse Weibull family
#' 
#' @author Amylkar Urrea Montoya, \email{amylkar.urrea@@udea.edu.co}
#' 
#' @description 
#' The Marshall-Olkin Extended Inverse Weibull family
#' 
#' @param mu.link defines the mu.link, with "log" link as the default for the mu parameter.
#' @param sigma.link defines the sigma.link, with "log" link as the default for the sigma.
#' @param nu.link defines the nu.link, with "log" link as the default for the nu parameter.
#' 
#' @seealso \link{dMOEIW}
#' 
#' @details 
#' The Marshall-Olkin Extended Inverse Weibull distribution with parameters \code{mu}, 
#' \code{sigma} and \code{nu} has density given by
#' 
#' \eqn{f(x) = \frac{\mu \sigma \nu x^{-(\sigma + 1)} exp\{{-\mu x^{-\sigma}}\}}{\{\nu -(\nu-1) exp\{{-\mu x ^{-\sigma}}\} \}^{2}},}
#' 
#' for x > 0. 
#' 
#' @returns Returns a gamlss.family object which can be used to fit a MOEIW distribution in the \code{gamlss()} function.
#' 
#' @example examples/examples_MOEIW.R 
#' 
#' @references
#' \insertRef{okasha2017}{RelDists}
#'
#' @importFrom Rdpack reprompt
#' @importFrom gamlss.dist checklink
#' @importFrom gamlss rqres.plot
#' @export
MOEIW <- function (mu.link="log", sigma.link="log", nu.link="log") {
  mstats <- checklink("mu.link", "Marshall-Olkin Extended Inverse Weibull", 
                      substitute(mu.link), c("log", "own"))
  dstats <- checklink("sigma.link", "Marshall-Olkin Extended Inverse Weibull",
                      substitute(sigma.link), c("log", "own"))
  vstats <- checklink("nu.link", "Marshall-Olkin Extended Inverse Weibull", 
                      substitute(nu.link), c("log", "own"))
  
  structure(list(family=c("MOEIW", "Marshall-Olkin Extended Inverse Weibull"), 
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
                   A    <- (nu - 1) * exp(-mu * y^(-sigma)) * y^(-sigma)
                   B    <- nu - (nu - 1) * exp(-mu * y^(-sigma))
                   dldm <- 1 / mu - y^(-sigma) - 2 * A / B
                   dldm
                 },
                 
                 dldd = function(y, mu, sigma, nu) {
                   C    <- 2 * (nu - 1) * exp(-mu * y^(-sigma)) * mu * log(y) * y^(-sigma)  
                   D    <- nu - (nu - 1) * exp(-mu * y^(-sigma)) 
                   dldd <- 1 / sigma - log(y) + mu * log(y) * y ^(-sigma) + C / D
                   dldd
                 },
                 
                 dldv = function(y, mu, sigma, nu) {
                   G    <- 2 * (1 - exp(-mu * y ^(-sigma)))
                   H    <- nu - (nu - 1) * exp(-mu * y^(-sigma))
                   dldv <- 1 / nu - G / H
                   dldv
                 },
                 
                 
                 # Segundas derivadas ---------------------------------
                 d2ldm2 = function(y, mu, sigma, nu) {
                   A    <- (nu - 1) * exp(-mu * y^(-sigma)) * y^(-sigma)
                   B    <- nu - (nu - 1) * exp(-mu * y^(-sigma))
                   dldm <- 1 / mu - y^(-sigma) - 2 * A / B
                   d2ldm2 <- -dldm * dldm
                   d2ldm2
                 },
                 
                 d2ldmdd = function(y, mu, sigma, nu) {
                   A    <- (nu - 1) * exp(-mu * y^(-sigma)) * y^(-sigma)
                   B    <- nu - (nu - 1) * exp(-mu * y^(-sigma))
                   dldm <- 1 / mu - y^(-sigma) - 2 * A / B
                   C    <- 2 * (nu - 1) * exp(-mu * y^(-sigma)) * mu * log(y) * y^(-sigma)  
                   D    <- nu - (nu - 1) * exp(-mu * y^(-sigma)) 
                   dldd <- 1 / sigma - log(y) + mu * log(y) * y ^(-sigma) + C / D
                   d2ldmdd <- -dldm * dldd
                   d2ldmdd
                 },
                 
                 d2ldmdv = function(y, mu, sigma, nu) {
                   A    <- (nu - 1) * exp(-mu * y^(-sigma)) * y^(-sigma)
                   B    <- nu - (nu - 1) * exp(-mu * y^(-sigma))
                   dldm <- 1 / mu - y^(-sigma) - 2 * A / B
                   G    <- 2 * (1 - exp(-mu * y ^(-sigma)))
                   H    <- nu - (nu - 1) * exp(-mu * y^(-sigma))
                   dldv <- 1 / nu - G / H
                   d2ldmdv <- -dldm * dldv
                   d2ldmdv
                 },
                 
                 d2ldd2 = function(y, mu, sigma, nu) {
                   C    <- 2 * (nu - 1) * exp(-mu * y^(-sigma)) * mu * log(y) * y^(-sigma)  
                   D    <- nu - (nu - 1) * exp(-mu * y^(-sigma)) 
                   dldd <- 1 / sigma - log(y) + mu * log(y) * y ^(-sigma) + C / D
                   d2ldd2 <- -dldd * dldd
                   d2ldd2
                 },
                 
                 d2ldddv = function(y, mu, sigma, nu) {
                   C    <- 2 * (nu - 1) * exp(-mu * y^(-sigma)) * mu * log(y) * y^(-sigma)  
                   D    <- nu - (nu - 1) * exp(-mu * y^(-sigma)) 
                   dldd <- 1 / sigma - log(y) + mu * log(y) * y ^(-sigma) + C / D
                   G    <- 2 * (1 - exp(-mu * y ^(-sigma)))
                   H    <- nu - (nu - 1) * exp(-mu * y^(-sigma))
                   dldv <- 1 / nu - G / H
                   d2ldddv <- -dldd * dldv
                   d2ldddv
                 },
                 
                 d2ldv2 = function(y, mu, sigma, nu) {
                   G    <- 2 * (1 - exp(-mu * y ^(-sigma)))
                   H    <- nu - (nu - 1) * exp(-mu * y^(-sigma))
                   dldv <- 1 / nu - G / H
                   d2ldv2 <- -dldv * dldv
                   d2ldv2
                 },
                 
                 G.dev.incr = function(y, mu, sigma, nu, ...) -2*dMOEIW(y, mu, sigma, nu, log=TRUE), 
                 rqres = expression(rqres(pfun="pMOEIW", type="Continuous", y=y, mu=mu, sigma=sigma, nu=nu)), 
                 
                 mu.initial = expression(mu       <- rep(1, length(y))), 
                 sigma.initial = expression(sigma <- rep(1, length(y))), 
                 nu.initial = expression(nu       <- rep(1, length(y))),
                 
                 mu.valid = function(mu)       all(mu > 0), 
                 sigma.valid = function(sigma) all(sigma > 0), 
                 nu.valid = function(nu)       all(nu > 0),
                 
                 y.valid = function(y) all(y >= 0)
  ), 
  class=c("gamlss.family", "family"))
}
